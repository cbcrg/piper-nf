/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'Piper-NF'.
 *
 *   Piper-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Piper-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Piper-NF.  If not, see <http://www.gnu.org/licenses/>.
 */


import com.google.common.collect.Multiset
import com.google.common.collect.HashMultiset
import nextflow.util.CacheHelper
import java.nio.file.Files
import java.util.regex.Pattern


/* 
 * Main Piper-NF pipeline script
 *
 * @authors
 * Giovanni Bussotti <giovannibussotti@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Pablo Prieto <poena.funesta@gmail.com> 
 *
 *
 * Defines the pipeline parameters.
 * The values in the 'params' map can be overridden on the command line by specifying a
 * option prefixed with a double '-' char, for example
 *
 * $ nextflow piper.nf --query=<path to your file name>
 *
 */

params.queryChunkSize = 1000
params.query = "${baseDir}/tutorial/5_RNA_queries.fa"
params.genomesDb = 'db'
params.resultDir = 'result'
params.blastStrategy = 'ncbi-blast'     // the blast tool to be used, choose between: ncbi-blast, wu-blast
params.alignStrategy = 'slow_pair'      // defines the T-Coffee alignment method
params.exonerateSuccess = '1'
params.exonerateMode = 1000             // max lines read by exonerate. Use 'exhaustive' to read all 
params.exonerateChunkSize = 2500
params.repeatCov = 20
params.extension = 20000
params.extensionFile = 'no'
params.osname = 'linux'


// these parameters are mutually exclusive
// Input genome can be specified by
// - genomes-file: a file containing the list of genomes FASTA to be processed
// - genomes-list: a comma separated list of genomes FASTA file
// - genomes-folder: a directory containing a folder for each genome FASTA file
params['genomes-file'] = null
params['genomes-list'] = null
params['genomes-folder'] = "${baseDir}/tutorial/genomes/"

queryFile = file(params.query)
dbPath = file(params.genomesDb)

if( !dbPath.exists() && !dbPath.mkdirs() ) {
    exit 1, "Cannot create genomes-db path: $dbPath -- check file system permissions"
}
log.info "Using genomes-db path: $dbPath"


/*
 * Verify that user has specified a valid BLAST strategy parameter
 */
assert params.blastStrategy in ['ncbi-blast','wu-blast']

/*
 * dump some info
 */

log.info "P I P E R - RNA mapping pipeline - ver 1.4.2"
log.info "============================================"
log.info "query               : ${queryFile}"
log.info "genomes-db          : ${dbPath}"
log.info "query-chunk-size    : ${params.queryChunkSize}"
log.info "genomes-file        : ${params.genomesFile ?: '-'}"
log.info "genomes-list        : ${params.genomesList ?: '-'}"
log.info "genomes-folder      : ${params.genomesFolder ?: '-'}"
log.info "result-dir          : ${params.resultDir}"
log.info "blast-strategy      : ${params.blastStrategy}"
log.info "align-strategy      : ${params.alignStrategy}"
log.info "exonerate-success:  : ${params.exonerateSuccess}"
log.info "exonerate-mode:     : ${params.exonerateMode}"
log.info "exonerate-chunk-size: ${params.exonerateChunkSize}"
log.info "repeat-cov          : ${params.repeatCov}"
log.info "\n"

/*
 * Find out all the genomes files in the specified directory.
 *
 * More in detail teh 'sourceGenomesPath' points to a directory having a
 * sub-folder for each genome it is required to process.
 *
 * Each sub-folder must contain the genome FASTA file to be processed.
 *
 * The sub-folder name is used to identify the genome in the computation.
 *
 * All the genomes names found in this path are put in a list named 'formatName',
 * which control the pipeline execution.
 *
 */

allGenomes = getGenomesMap(params)

/*
 * Split the query input file in many small files (chunks).
 *
 * The number of sequences in each chunk is controlled by the parameter 'queryChunkSize'
 * The chunk files are saved in a local folder define by the variable 'querySplits'
 *
 */

// create a folder that may be cached, using the 'queryFile' and the number chunks as cache key
querySplits = cacheableDir([queryFile, params.queryChunkSize])

if( querySplits.empty() ) {
    log.info "Splitting query file: $queryFile .."
    chunkCount=0
    queryFile.splitFasta( by: params.queryChunkSize ) {
        querySplits.resolve("seq_${chunkCount++}").text = it
    }
    log.info "Created $chunkCount input chunks to path: ${querySplits}"
}
else {
    log.info "Cached query splits > ${querySplits.list().size()} input query chunks"
}


// use a set since there should be not repetition
allQueryIDs = new HashSet()

queryEntries = cacheableDir(queryFile)
log.debug "Queries entries path: $queryEntries"

queryFile.splitFasta(record: [header:true, text:true]) { record ->
    String queryId = record.header.replaceAll( /^([^\s\|]*).*$/, '$1' )
    log.debug "Query entry id: $queryId"

    allQueryIDs << queryId
    // store the chunk to a file named as the 'queryId'
    def fileEntry = queryEntries.resolve(queryId)
    if( fileEntry.isEmpty() ) {
        fileEntry.text = record.text
    }
}

/*
 * Create the required databases (BLAST,CHR) if they do not exists.
 *
 * This task is executed for each genome in the list 'formatName'
 * The tasks 'sends' out the name of the genome to be processed
 * by the next step in the pipeline using the variable 'blastName'
 *
 */
 
def sed_cmd = (params.osname == 'osx' ? 'gsed' : 'sed')
def split_cmd = (params.osname == 'osx' ? 'gcsplit' : 'csplit') 

/*
 * Creates two channels that feed the 'formatBlast' and 'formatChr' processes.
 * Both of there emits tuples like (species, genome fasta file)
 */

Channel
        .from(allGenomes.entrySet())
        .map { entry -> [ entry.key, entry.value ] }
        .into{ fmtChrParams; fmtBlastParams }

/*
 * Given the genome FASTA file, format it to the BLAST binary format
 */
process formatBlast {

    storeDir { dbPath / species }

    input:
    set (species, file(genome_fa)) from fmtBlastParams

    output:
    set (species, file { blast_db }) into fmtBlastOut

    script:
    blast_db = "${params.blastStrategy}-db"

    if( params.blastStrategy == 'ncbi-blast' )
    """
    mkdir -p $blast_db
    makeblastdb -dbtype nucl -in ${genome_fa} -out $blast_db/db
    """

    else if( params.blastStrategy == 'wu-blast' )
    """
    mkdir -p $blast_db
    xdformat -n -o ${blast_db}/db ${genome_fa}
    """
}

/*
 * Given the genome FASTA file splits in a file for each sequence
 */
process formatChr {

    storeDir { dbPath / species }

    input:
    set (species, file(genome_fa)) from fmtChrParams

    output:
    set (species, file('chr')) into fmtChrOut

    script:
    """
    ## split the fasta in a file for each sequence 'seq_*'
    ${split_cmd} ${genome_fa} '%^>%' '/^>/' '{*}' -f seq_ -n 5

    ## create the target folder
    mkdir -p chr

    ## rename and move to the target folder
    for x in seq_*; do
    SEQID=`grep -E "^>" \$x | ${sed_cmd} -r 's/^>(\\S*).*/\\1/'`
    mv \$x chr/\$SEQID;
    done
    """
}


/*
 * Implements the BLAST step
 */

fmtBlastOut.spread( querySplits.listFiles() ).set { blast_in }

process blast {
    input:
    set (blastId, file(blast_db), file('blastQuery') ) from blast_in

    output:
    set (blastId, 'blastQuery', '*.mf2' ) into blast_result

    script:

    if( params.blastStrategy == 'ncbi-blast' )

        """
        fmt='6 qseqid sseqid evalue score qgi bitscore length nident positive mismatch pident ppos qacc gaps gaopen qaccver qlen qframe qstart qend sframe sstart send'
        blastn -db $blast_db/db -query blastQuery -outfmt "\$fmt" -num_threads ${task.cpus} > ${blastId}.mf2
        """


    else if( params.blastStrategy == 'wu-blast' )
        """
        wu-blastn $blast_db/db blastQuery -mformat=2 -e 0.00001 -cpus ${task.cpus} -filter=seg -lcfilter -errors -novalidctxok -nonnegok > ${blastId}.mf2
        """

}

/*
 * Split the blast result in chunk of at most *params.exonerateChunkSize* lines
 */

blast_result
    .splitText( by: params.exonerateChunkSize, elem: 2 )
    .set { blast_chunks }

/*
 * Join together the output of 'formatChr' step with the 'blast' step
 *
 * Channel 'fmtChrOut' emits tuples with these elements (species, chr_db)
 * Channel 'blast_result' emits tuples with these elements ( species, blast_query, blast_hits )
 *
 * Finally, the channel 'exonerate_in' emits ( species, chr_db, blast_query, blast_hits )
 */

fmtChrOut
     .cross( blast_chunks )
     .map { chr, blast ->
            [ chr[0], chr[1], blast[1], blast[2] ]
     }
     .set { exonerate_in }


/*
 * Collect the BLAST output chunks and apply the 'exonerate' function
 */

process exonerate {
    input:
    set ( species, file(chr_db), file(exonerateQuery), file('blast_result') ) from exonerate_in

    output:
    set ( species, '*.fa', '*.gtf') into exonerate_out

    """
    exonerateRemapping.pl \
        -query ${exonerateQuery} \
        -mf2 blast_result \
        -targetGenomeFolder ${chr_db} \
        -exonerate_lines_mode ${params.exonerateMode} \
        -exonerate_success_mode ${params.exonerateSuccess} \
        -ner no \
	-extensionFile ${params.extensionFile} \
	-extension ${params.extension}

    if [ -s blast_result.fa ]; then
      repeat.pl blast_result.fa blast_result.ex.gtf ${params.repeatCov}
      if [ -s rep${params.repeatCov}.fa ]; then
          mv blast_result.fa chunk.seq
          mv blast_result.ex.gtf chunk.ex.annot
          mv rep${params.repeatCov}.fa blast_result.fa
          mv rep${params.repeatCov}.ex.gtf blast_result.ex.gtf
      fi
    fi

    mv blast_result.fa '${species}.fa'
    mv blast_result.ex.gtf '${species}.ex.gtf'
    """
}

/*
 * Post-process 'exonerate' result
 *
 * It collects all the fasta and gtf files produced by the 'exonerate' step and the 'hitxx' suffix
 * in such a way that the 'xx' number is unique across the same species and queryId
 *
 */

Multiset hitSet = HashMultiset.create()

process normExonerate {

  input:
  set ( species, fasta, gtf ) from exonerate_out

  output:
  val normalizedGtf
  val normalizedFasta

  exec:
    def chunks = []
    def replace = []

    fasta.splitFasta(record:[id:true, sequence:true]) { record ->

        // parse the sequence id
        def matcher = (record.id =~ /(.*)_(hit\d*)(.*)/ )
        def (queryId, hitName, extra) = matcher[0][1..3]

        // create a multi-fasta file for each 'queryId'
        if( !allQueryIDs.contains(queryId) ) {
            println "Skipping queryId: $queryId -- since it's not contained in the source query"
            return
        }

        log.debug "normExonerate > Processing queryId: ${queryId}"

        // update the hit name
        def key = [species, queryId]
        def count = hitSet.add(key, 1) +1
        def newHit = "hit$count"
        if( hitName != newHit ) {
            log.debug "Replacing hitName: $hitName with: $newHit using key: $key"
            replace << [queryId: queryId, oldHit: hitName, newHit: newHit ]
            hitName = newHit
        }

        // now append the query content
        def item = ''
        item += ">${queryId}_${hitName}${extra}_${species}\n"
        item += record.sequence

        // return a pair made up of the queryId and the
        chunks << [ queryId, item ]
    }

    // clone the chunks list and assign it as an output value
    normalizedFasta = chunks.clone()

    if( !replace ) {
        normalizedGtf = gtf
        return
    }

    // normalizing hitNames
    def str = gtf.text
    replace.each {
        log.debug "normExonerate > Replacing hitName: $it in GTF file: $gtf"
        def pattern = Pattern.quote("hitName \"${it.queryId}_${it.oldHit}\";")
        str = str.replaceAll( ~/$pattern/, "hitName \"${it.queryId}_${it.newHit}\";" )
    }

    // creates a new gtf with the same name (in a different directory)
    normalizedGtf = workDir.resolve(gtf.name)
    normalizedGtf.text = str

}

/*
 *  Groups together all the sequences (in the fasta file) output by exonerate that has the same 'queryId'
 *  and save them in a file having the name $queryId.mfa
 *
 *  All of them are stores in the cacheable path 'fastaCacheDir'
 *
 */

queryFiles = [:]
queryEntries.eachFile { queryFiles["${it.name}.mfa"] = it }

normalizedFasta
        .flatMap()
        .collectFile( seed: queryFiles ) { queryId, sequence -> ["${queryId}.mfa", sequence] }
        .set { alignFasta }

/*
 * Aligns the the query sequence with hit sequences returned by exonerate
 */

process similarity {

    input:
    file seq from alignFasta

    output:
    file '*.sim' into similarity

    """
    t_coffee -in '${seq}' -method ${params.alignStrategy} -n_core ${task.cpus}
    t_coffee -other_pg seq_reformat -in *.aln -output sim > '${seq.baseName}.sim'
    """
}


similarity
    .collectFile { [it.baseName, it] }
    .toSortedList()
    .set { similarityFiles }



/*
 * Copy the GFT files produces by the Exonerate steps into the result (current) folder
 */

resultDir = file(params.resultDir)
resultDir.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
}

normalizedGtf.collectFile(storeDir: resultDir)

/*
 * Compute the similarity Matrix
 */
process matrix {
    echo true

    input:
    file queryFile
    file dbPath
    file similarityFiles

    output:
    file simMatrix

    """
    echo '\n====== Pipe-R sim matrix ======='
    mkdir data
    mv ${similarityFiles} data
    sim2matrix.pl -query $queryFile -data_dir data -genomes_dir $dbPath | tee simMatrix
    echo '\n'
    """
}

/*
 * Save the similarity matrix
 */
simMatrix.subscribe { file ->
    def target = resultDir.resolve('simMatrix.csv')
    println "\nSaving similarity matrix to file:\n> $target"
    file.copyTo(target)
}


// ----==== utility methods ====----


def parseGenomesFile( sourcePath ) {

    def result = [:]

    // parse the genomes input file files (genome-id, path to genome file)
    int count=0
    sourcePath.eachLine { line ->

        def genomeId
        def path

        def items = line.trim().split(/\s+/)
        if( items.size() > 1 ) {
            count++
            (path, genomeId) = items
        }
        else if( items.size() ==1 && items[0] ){
            count++
            genomeId = "gen${count}"
            path = items[0]
        }
        else {
            return
        }

        result[ genomeId ] = file(path)
    }

    result
}



def parseGenomesList(String genomesList) {

    def count=0
    def files = genomesList.split(',').collect { file(it.trim()) }
    def result = [:]

    files.each { genomeFile ->
        def genomeId = "gen${++count}"
        result[ genomeId ] = genomeFile
    }
    result
}

def parseGenomesFolder(sourcePath) {

    def result = [:]

    sourcePath.eachDir { path ->
        def fasta = path.listFiles().find { file -> file.name.endsWith('.fa') }
        if( fasta ) {
            result[ path.name ] = fasta
        }
    }
    result
}

def getGenomesMap(Map params) {

    def allGenomes = [:]
    // when the provided source path is a FILE
    // each line represent the path to a genome file
    if( params['genomes-file'] ) {
        def genomesFile = file(params['genomes-file'])
        if( genomesFile.empty() ) {
            exit 1, "Not a valid input genomes descriptor file: ${genomesFile}"
        }

        allGenomes = parseGenomesFile(genomesFile)
    }

    else if( params['genomes-list'] ) {
        allGenomes = parseGenomesList(params['genomes-list'])
    }

    else if( params['genomes-folder'] ) {
        def sourcePath = file(params['genomes-folder'])
        if( !sourcePath.exists() || sourcePath.empty() ) {
            exit 4, "Not a valid input genomes folder: ${sourcePath}"
        }

        allGenomes = parseGenomesFolder(sourcePath)
    }

    else {
        exit 5, "No input genome(s) provided -- Use one of the following CLI options 'genomes-file' or 'genomes-list' or 'genomes-folder' "
    }

    if( !allGenomes ) {
        exit 6, "No genomes found in path"
    }

    allGenomes.each { name, genome_fa ->
        log.info "Validating genome: $name -- file: ${genome_fa}"
        if( !genome_fa.exists() ) {
            exit 3, "Missing genome file: ${genome_fa}"
        }
    }

    return allGenomes
}

// ----===== TEST ====-------

def void testParseGenomesFile() {

    def source = file('test-source')
    try {
        source.text =
            '''
            x/file1.fa
            y/file2.fa   genx
            z/file3.fa
            '''

        def result = parseGenomesFile(source)

        assert result.size() == 3

        assert result['gen1'] == file('x/file1.fa')
        assert result['genx'] == file('y/file2.fa')
        assert result['gen3'] == file('z/file3.fa')

    }
    finally {
        source.delete()
    }
}



def void testParseGenomesList() {

      // call the function to test
      def result = parseGenomesList('alpha.fa, beta.fa, delta.fa')

      // verify result
      assert result.size() == 3

      assert result['gen1'] == file('alpha.fa')
      assert result['gen2'] == file('beta.fa')
      assert result['gen3'] == file('delta.fa')

}

def void testParseGenomesFolder() {

  def root = file('test-folder')

  try {
      // create the structure to test
      def folder1 = root.resolve('alpha')
      def folder2 = root.resolve('beta')
      def folder3 = root.resolve('delta')
      folder1.mkdirs()
      folder2.mkdirs()
      folder3.mkdirs()

      folder1.resolve('gen1.fa').text = 'uno'
      folder2.resolve('gen2.fa').text = 'due'
      folder3.resolve('gen3.fa').text = 'tre'

      // call the function to test
      def result = parseGenomesFolder(root)

      // verify result
      assert result.size() == 3

      assert result['alpha'] == folder1.resolve('gen1.fa')
      assert result['beta'] == folder2.resolve('gen2.fa')
      assert result['delta'] == folder3.resolve('gen3.fa')

  }
  finally {
     root.deleteDir()
  }


}
