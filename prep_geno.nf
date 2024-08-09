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

params.geno_dir = null
params.out_dir = "result"
params.max_maf = 0.9999
params.min_maf = 0.0001
params.miss = 0.05
params.hwe = 0.00001
params.keep = null
params.regions = null
params.chr = null
params.help = false


// Print usage and help

if (params.help) {
    log.info ''
    log.info 'D I G E N I C - N F'
    log.info '======================================================================='
    log.info 'Pre-process and filter UKB files, annotate them and generate pairs'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow run prep_geno.nf [options]'
    log.info ''
    log.info 'Parameters:'
    log.info " --geno_dir GENOTYPE_DIR     directory containing genotypes in BGEN format (default: $params.geno_dir)"
    log.info " --out_dir OUTPUT_DIR        output directory (default: $params.out_dir)"
    log.info " --min_maf MIN_MAF           minimum MAF allowed (default: $params.min_maf)"
    log.info " --max_maf MAX_MAF           maximum MAF allowed (default: $params.max_maf)"
    log.info " --miss MAX_MISS             maximum missingness (at variant level) allowed (default: $params.miss)"
    log.info " --hwe HWE                   hardy-weinberg P-value threshold (default: $params.hwe)"
    log.info " --keep KEEP                 single-column file with the list of individuals to keep (default: $params.keep)"
    log.info " --regions REGIONS           BED file with the regions to keep (default: $params.regions)"
    log.info " --chr CHR                   chr (single chr, comma-separated list, range via colon) (default: $params.chr)"
    log.info ''
    exit(1)
}

// Print parameter selection

log.info ''
log.info 'Parameters'
log.info '------------------'
log.info "Genotypes                    : ${params.geno_dir}"
log.info "Output directory             : ${params.out_dir}"
log.info "Minimum MAF                  : ${params.min_maf}"
log.info "Maximum MAF                  : ${params.max_maf}"
log.info "Max. missingess rate (vars)  : ${params.miss}"
log.info "HWE P-value treshold         : ${params.hwe}"
log.info "List of individuals          : ${params.keep}"
log.info "List of regions (BED)        : ${params.regions}"
log.info "List of chromosomes          : ${params.chr}"
log.info ''


// Mandatory options

if (!params.geno_dir) {
    exit 1, 'Genotype directory not specified.'
} else if (! params.keep) {
    exit 1, 'Individual list file not specified'
} else if (! params.regions) {
    exit 1, 'Regions BED file not specified'
} else if (! params.chr) {
    exit 1, 'Chromosome list not specified'
} 

// Expand chr parameter

chrlist = []
if (params.chr =~ /,/) {
    chrlist = params.chr.tokenize(',')
} else {
    chrlist = params.chr
}

// Processes

process Filter {

    tag { "chr${chr}" }

    input:
    tuple val(chr), file(bgen), file(sample)
    file(keep)
    file(regions)
 
    output:
    tuple file("chr${chr}.pgen"), file("chr${chr}.pvar"), file("chr${chr}.psam")

    """
    plink2 --bgen $bgen ref-first --sample $sample \
           --maf $params.min_maf --max-maf $params.max_maf \
           --geno $params.miss \
           --hwe $params.hwe midp \
           --keep $params.keep \
           --extract bed1 $params.regions \
           --remove-nosex \
           --rm-dup exclude-all \
           --make-pgen --out chr$chr --threads 12
    """
}

process Merge {

   publishDir "${params.out_dir}", mode: 'copy'

   input:
   file(f)
 
   output:
   tuple file('aid.pgen'), file('aid.psam'), file('aid.pvar')   

   """
   ls chr*pgen | sort -V | sed 's/.pgen//' > tomerge.txt
   plink2 --pmerge-list tomerge.txt --out aid
   """
}

process Pairs {

   publishDir "${params.out_dir}", mode: 'copy'

   input:
   tuple file(pgen), file(psam), file(pvar)


   output:
   file('pairs.txt')

   """
   awk 'NR>1{print \$1":"\$2"-"\$2}' $pvar > variants.txt
   combos_diffchr.py -i variants.txt > pairs.txt
   """
}

// Pipeline

workflow {
In = Channel.of(chrlist)
        .flatten()
        .map { it -> [it, file("${params.geno_dir}/ukb22828_c${it}_b0_v3.bgen"), file("${params.geno_dir}/ukb22828_c${it}_b0_v3_*.sample")]}
aid_plink2 = Filter(In, file(params.keep), file(params.regions)) | flatten | collectFile | collect | Merge
Pairs(aid_plink2)
}
