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

params.queryChunkSize = 100
params.query = 'tutorial/5_RNA_queries.fa'
params.genomesDb = 'db'
params.resultDir = 'result'
params.blastStrategy = 'ncbi-blast'     // the blast tool to be used, choose between: ncbi-blast, wu-blast
params.alignStrategy = 'slow_pair'      // defines the T-Coffee alignment method
params.exonerateSuccess = '1'
params.exonerateMode = 'exhaustive'


// these parameters are mutually exclusive
// Input genome can be specified by
// - genomes-file: a file containing the list of genomes FASTA to be processed
// - genomes-list: a comma separated list of genomes FASTA file
// - genomes-folder: a directory containing a folder for each genome FASTA file
params['genomes-file'] = null
params['genomes-list'] = null
params['genomes-folder'] = "tutorial/genomes/"

queryFile = file(params.query)
dbPath = file(params.genomesDb)

if( !dbPath.exists() ) {
    log.warn "Creating genomes-db path: $dbPath"
    if( !dbPath.mkdirs() ) {
        exit 1, "Cannot create genomes-db path: $dbPath -- check file system permissions"
    }
}

log.info "P I P E R - RNA mapping pipeline"
log.info "================================"
log.info "query              : ${queryFile}"
log.info "genomes-db         : ${dbPath}"
log.info "query-chunk-size   : ${params.queryChunkSize}"
log.info "result-dir         : ${params.resultDir}"
log.info "blast-strategy     : ${params.blastStrategy}"
log.info "align-strategy     : ${params.alignStrategy}"
log.info "exonerate-success  : ${params.exonerateSuccess}"
log.info "exonerate-mode     : ${params.exonerateMode}"
log.info "pool-size          : ${config.poolSize}"
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

allGenomes = [:]

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
    queryFile.chunkFasta( params.queryChunkSize ) { sequence ->
        def file = querySplits.resolve("seq_${chunkCount++}")
        file.text = sequence
    }
    log.info "Created $chunkCount input chunks to path: ${querySplits}"
}
else {
    log.info "Cached query splits > ${querySplits.list().size()} input query chunks"
}


allQueryIDs = []
queryFile.chunkFasta() { String chunk ->
    allQueryIDs << chunk.readLines()[0].replaceAll( /^>(\S*).*$/, '$1' )
}

/*
 * Create the required databases (BLAST,CHR) if they does not exists.
 *
 * This task is executed for each genome in the list 'formatName'
 * The tasks 'sends' out the name of the genome to be processed
 * by the next step in the pipeline using the variable 'blastName'
 */


def sed_cmd = (System.properties['os.name'] == 'Mac OS X' ? 'gsed' : 'sed')
def split_cmd = (System.properties['os.name'] == 'Mac OS X' ? 'gcsplit' : 'csplit')

/*
 * creates two channels that feed the 'formatBlast' and 'formatChr' processes
 */
fmtBlastParams = Channel.create()
fmtChrParams = Channel.create()
Channel
        .from(allGenomes.entrySet())
        .map { entry -> [ entry.key, entry.value ] }   // <-- return a tuple like (specie, genome_fasta)
        .split( fmtChrParams, fmtBlastParams )

/*
 * Given the genome FASTA file, format to a binary format using the BLAST
 */
process formatBlast {

    storeDir dbPath

	input:
    set (specie, genome_fa) from fmtBlastParams

    output:
    set (specie, file("$blast_db")) into fmtBlastOut

    script:
    blast_db = "$specie/${params.blastStrategy}-db"
    """
    ## Create the target folder
    mkdir -p $blast_db

    ## Format the BLAST DB
    x-format.sh ${params.blastStrategy} ${genome_fa} ${blast_db}
    """
}

/*
 * Given the genome FASTA file splits in a file file for each sequence
 */
process formatChr {

    storeDir dbPath

	input:
    set (specie, genome_fa ) from fmtChrParams

    output:
    set (specie, file("$chr_db")) into fmtChrOut

    script:
    chr_db = "${specie}/chr"
    """
    ## split the fasta in a file for each sequence 'seq_*'
    ${split_cmd} ${genome_fa} '%^>%' '/^>/' '{*}' -f seq_ -n 5

    ## create the target folder
    mkdir -p ${chr_db}

    ## rename and move to the target folder
    for x in seq_*; do
    SEQID=`grep -E "^>" \$x | ${sed_cmd} -r 's/^>(\\S*).*/\\1/'`
    mv \$x ${chr_db}/\$SEQID;
    done
    """
}


/*
 * Implements the BLAST step
 */

blast_in = fmtBlastOut.spread( querySplits.listFiles() )

process blast {
    input:
    set (blastId, file(blast_db), file('blastQuery') ) from blast_in

    output:
    set (blastId, 'blastQuery', '*.mf2' ) into blast_result

    script:

    if( params.blastStrategy == 'ncbi-blast' )

        """
    	fmt='6 qseqid sseqid evalue score qgi bitscore length nident positive mismatch pident ppos qacc gaps gaopen qaccver qlen qframe qstart qend sframe sstart send'
    	blastn -db $blast_db/db -query blastQuery -outfmt "\$fmt" > ${blastId}.mf2
    	"""


    else if( params.blastStrategy == 'wu-blast' )
        """
    	wu-blastn $dblast_db/db blastQuery -mformat=2 -e 0.00001 -cpus 1 -filter=seg -lcfilter > ${blastId}.mf2
    	"""

}


exonerate_in = fmtChrOut                           // fmtChrOut: emits  (specie, chr_db)
                 .cross( blast_result )           // blast_result: emits (specie, blast_query, blast_hits )
                 .map { chr, blast ->
                        [ chr[0], chr[1], blast[1], blast[2] ]   // returns ( specie, chr_db, blast_query, blast_hits )
                      }


/*
 * Collect the BLAST output chunks and apply the 'exonerate' function
 */

process exonerate {
    input:
    set ( specie, file(chr_db), file(exonerateQuery), file(blastResult) ) from exonerate_in

    output:
    file '*.fa' into exonerateOut mode flatten
    file '*.gtf' into exonerateGtf mode flatten

    """
    ## apply exonerate
    exonerateRemapping.pl -query ${exonerateQuery} -mf2 $blastResult -targetGenomeFolder $chr_db -exonerate_lines_mode ${params.exonerateMode} -exonerate_success_mode ${params.exonerateMode} -ner no
    if [ ! -s *.fa ]; then exit 0; fi

    ## exonerateRemapping create a file named '*.fa'
    ## split the exonerate result into single files
    ${split_cmd} *.fa '%^>%' '/^>/' '{*}' -f .seq_ -n 5
    mv *.fa .exonerate.fa

    ## rename the seq_xxx files so that the file name match the seq fasta id
    ## plus append the specie to th sequence id
    for x in .seq_*; do
      SEQID=`grep '>' \$x`
      FILENAME=`grep '>' \$x | ${sed_cmd} -r 's/^>(.*)_hit\\d*.*\$/\\1/'`
      printf "\${SEQID}_${specie}\\n" >> \${FILENAME}.fa
      cat \$x | grep -v '>' >> \${FILENAME}.fa
    done

    """
}



fastaToMerge = exonerateOut.filter { file -> file.baseName in allQueryIDs  }

process prepare_mfa {
    merge true

    input:
    file fastaToMerge

    output:
    file '*.mfa' into fastaToAlign mode flatten

    """
    # Extract the file name w/o the extension
    baseName="${fastaToMerge.baseName}"

    # Only the first time append the query sequence
    if [ ! -e \$baseName.mfa ]; then
    perl -n -e '\$on=(/^>('\$baseName')\$/) if (/^>/); print \$_ if (\$on);' $queryFile > \$baseName.mfa
    fi

    # Append the exonerate result
    cat $fastaToMerge >> \$baseName.mfa
    """
}

process align {
    input:
    file fastaToAlign

    output:
    file '*.aln' into alignment

    """
    t_coffee -in $fastaToAlign -method ${params.alignStrategy} -n_core 1
    """
}

process similarity {
    merge true

    input:
    file alignment

    output:
    file '*' into similarity

    """
    t_coffee -other_pg seq_reformat -in $alignment -output sim > ${alignment.baseName}
    """
}

/*
 * Copy the GFT files produces by the Exonerate steps into the result (current) folder
 */

resultDir = file(params.resultDir)
resultDir.with {
    if( !empty() ) { deleteDir() }
    mkdirs()
}

exonerateGtf.each { file ->
    if( file.size() == 0 ) return
    def name = file.name
    def gtfFileName = resultDir.resolve(name)
    gtfFileName << file.text
}



/*
 * Compute the similarity Matrix
 */
process matrix {
    echo true

    input:
    file similarity

    output:
    file 'simMatrix'

    """
    echo '\n====== Pipe-R sim matrix ======='
    mkdir data
    mv ${similarity} data
    sim2matrix.pl -query $queryFile -data_dir data -genomes_dir $dbPath | tee simMatrix
    echo '\n'
    """
}

simMatrix.subscribe { file ->
    file.copyTo( resultDir.resolve('simMatrix.csv') )
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

      def db = file('db')

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
