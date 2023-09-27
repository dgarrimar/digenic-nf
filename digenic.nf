#!/bin/env nextflow

/*
 * Copyright (c) 2023, Diego Garrido-Martín
 *
 * Simulation setting to benchmark the TIE and power of multivariate 
 * methods (MANTA, MANOVA and GEMMA) in the context of multi-trait 
 * GWAS studies.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */


//  Define parameters

params.help = false

// Print usage and help

if (params.help) {
    log.info ''
    log.info 'D I G E N I C - N F'
    log.info '======================================================================='
    log.info 'Compute goodness-of-fit chi-square statistic for a set of variant pairs'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run digenic.nf [options]'
    log.info ''
    log.info 'Parameters:'
    log.info ' --geno GENOTYPES           genotypes in VCF format (default: $geno)'
    log.info ''
    exit(1)
}

// Print parameter selection

log.info ''
log.info 'Parameters'
log.info '------------------'
log.info "Genotypes                    : ${params.geno}"
log.info ''

// Mandatory options

if (!params.geno) {
    exit 1, "Genotype file not specified."
} 

// Pipeline

/*
 *  Generate pairs 
 */

process pairs {

    input:
    path(geno)

    output:
    path('pairs.txt')    

    script:
    """
    bcftools query -f '%CHROM:%POS\\n' $geno > variants.txt
    combos.py -i variants.txt > pairs.txt
    """
}

process chi2 {

   tag {id}

   input:
   tuple val(id), path(chunk)
   path(geno) 
   path(geno_index)

   output:
   tuple val(id), path('sstats.*')

   shell:
   ''' 
   one_bkp=''
   cat !{chunk} | while read line; do 
       one=$(echo $line | cut -d ' ' -f1)
       two=$(echo $line | cut -d ' ' -f2)
       if [[ $one == $one_bkp ]]; then
           bcftools view -r $two -H !{geno} | cut -f10- > two.txt
       else
           one_bkp=$one
           bcftools view -r $one -H !{geno} | cut -f10- > one.txt
           bcftools view -r $two -H !{geno} | cut -f10- > two.txt
       fi
       cat one.txt two.txt > vp.txt
       touch sstats.!{id}
       break
   done
   '''
}

workflow {
    geno = channel.value(params.geno) 
    geno_index = geno.map{it.replace(".vcf.gz", ".vcf.gz.tbi")}
    chunks = pairs(geno) | splitText(by: 500, file: 'chunk') | map {[it.name.replace("chunk.",""), it]}  
    chunks.view()
    chi2 (chunks, geno, geno_index) \
      | view
}
