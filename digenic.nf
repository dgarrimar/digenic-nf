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
    nchoose2.py -n \$(bcftools view -H $geno | wc -l) -o pairs.txt
    """
}

workflow {
    channel.fromPath(params.geno) \
      | pairs \
      | view
}
