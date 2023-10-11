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

params.dir = 'result'
params.out = 'sstats.txt'
params.geno = null
params.pairs = null
params.cs = 20
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
    log.info " --geno GENOTYPES            genotypes in VCF format (default: $params.geno)"
    log.info " --pairs PAIRS               variant pairs (default: $params.pairs)"
    log.info " --cs CHUNK SIZE             chunk size (default: $params.cs)"
    log.info " --dir DIRECTORY             output directory (default: $params.dir)"
    log.info " --out OUTPUT                output file (default: $params.out)"
    log.info ''
    exit(1)
}

// Print parameter selection

log.info ''
log.info 'Parameters'
log.info '------------------'
log.info "Genotypes                    : ${params.geno}"
log.info "Variant pairs                : ${params.pairs}"
log.info "Chunk size                   : ${params.cs}"
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
log.info ''


// Mandatory options

if (!params.geno) {
    exit 1, 'Genotype file not specified.'
} 

// Processes

process prep_vcf {

    input:
    path(geno)

    output:
    tuple path('gt.traw.gz'), path('gt.traw.gz.tbi')

    """
    plink2 --vcf $geno --recode Av --out gt
    cut -f3,5,6 --complement gt.traw > tmpfile; mv tmpfile gt.traw
    bgzip gt.traw
    tabix -s 1 -b 3 -e 3 -S 1 gt.traw.gz
    """
}

process prep_plink {

    input:
    tuple path(pgen),path(pvar),path(psam)

    output:
    tuple path('gt.traw.gz'), path('gt.traw.gz.tbi')

    """
    geno=\$(echo $pgen | sed -r 's/\\..+//')
    plink2 --pfile \$geno -make-bed --out \$geno
    plink2 --bfile \$geno --recode Av --out gt
    cut -f3,5,6 --complement gt.traw > tmpfile; mv tmpfile gt.traw
    bgzip gt.traw
    tabix -s 1 -b 3 -e 3 -S 1 gt.traw.gz
    """


}

process pairs {

    input:
    path(geno)

    output:
    path('pairs.txt')    

    script:
    """
    bcftools query -f '%CHROM:%POS-%POS\\n' $geno > variants.txt
    combos.py -i variants.txt > pairs.txt
    """
}

process chi2 {

   tag { id }

   input:
   tuple val(id), path(chunk)
   tuple path(genop), path(genop_idx)

   output:
   path('sstats.*')

   shell:
   ''' 
   one_bkp=''
   cat !{chunk} | while read line; do 
       one=$(echo $line | cut -d ' ' -f1)
       two=$(echo $line | cut -d ' ' -f2)
       if [[ $one == $one_bkp ]]; then
           tabix !{genop} $two | cut -f4- | sed 's/\\t/\\n/g' > two.txt
       else
           one_bkp=$one
           tabix !{genop} $one | cut -f4- | sed 's/\\t/\\n/g' > one.txt
           tabix !{genop} $two | cut -f4- | sed 's/\\t/\\n/g' > two.txt
       fi
       paste one.txt two.txt > all
       echo -e "$one $two $(chi2.R)"
   done > sstats.!{id}
   '''
}

process end {
    publishDir "${params.dir}", mode: 'copy'

    input:
    path(sstats)

    output:
    path(sstats)

    """
    """
}

// Pipeline

workflow {
    // geno = channel.value(params.geno) 
    geno = Channel.value([file("${params.geno}.pgen"), file("${params.geno}.pvar"), file("${params.geno}.psam") ])
    genop = prep_plink(geno)
    if( params.pairs ) {
        chunks = Channel.fromPath(params.pairs).splitText(by: params.cs, file: "chunk").map(){[it.name.replace("chunk.", ""), it]}
    } else {
        chunks = pairs(geno) | splitText(by: params.cs)
    }
    chi2 (chunks, genop) | collectFile(name: "${params.out}", sort: {it.name}) | end 
}
