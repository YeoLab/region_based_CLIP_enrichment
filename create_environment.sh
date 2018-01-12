#!/bin/bash
conda create -y --name regionnormalize python=2.7;

source activate regionnormalize;
pip install cwlref-runner;
conda install -y samtools;
conda install -c r r-essentials;
conda install -y perl;
conda install -y -c conda-forge perl-app-cpanminus;
cpanm install Statistics::Basic;
cpanm install Statistics::Distributions;
cpanm install Statistics::R;

export PATH=${PWD}/bin:$PATH
export PATH=${PWD}/bin/perl:$PATH
export PATH=${PWD}/cwl:$PATH
export PATH=${PWD}/wf:$PATH

