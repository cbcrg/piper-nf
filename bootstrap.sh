#!/bin/bash
#
#  Copyright (c) 2013, Centre for Genomic Regulation (CRG) and the authors 
#
#  This file is part of Piper-NF.
#
#  Piper-NF is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  Piper-NF is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Piper-NF.  If not, see <http://www.gnu.org/licenses/>.
#
#

set -e

install() {

  #
  # Install missing packages 
  # 
  sudo apt-get update --fix-missing
  sudo apt-get install -y openjdk-7-jre-headless wget
  sudo apt-get install -y ncbi-blast+


  #
  # T-COFFEE 
  # 
  wget -q -r -l1 --no-parent -A "T-COFFEE_installer*_linux_x64.bin" http://tcoffee.org/Packages/Stable/Latest/linux/
  tcoffee=`ls tcoffee.org/Packages/Stable/Latest/linux/*.bin`
  chmod +x $tcoffee
  $tcoffee --mode unattended --user_email tcoffee.msa@gmail.com
  rm -rf tcoffee.org

  # 
  # Exonerate 
  # 
  wget -q http://www.ebi.ac.uk/~guy/exonerate/exonerate-2.2.0-x86_64.tar.gz
  tar xf exonerate-2.2.0-x86_64.tar.gz 
  exonerate=`ls -d $HOME/exonerate-*/bin`
  printf "export PATH=$exonerate:\$PATH\n" >> ~/.profile
  rm -rf exonerate-2.2.0-x86_64.tar.gz
  
  #
  # chr_subseq
  #
  mkdir "$HOME/bin"
  wget -q -O $HOME/bin/chr_subseq http://tcoffee.org/Packages/Archive/chr_subseq
  chmod +x $HOME/bin/chr_subseq  
  printf "export PATH=\$HOME/bin:\$PATH\n" >> ~/.profile
  
} 

# Exit if already bootstrapped.
test -f /etc/bootstrapped && exit

export -f install
su vagrant -c 'install'

# Mark as bootstrapped 
date > /etc/bootstrapped