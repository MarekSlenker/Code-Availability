#!/bin/bash

# Receive nucleotide substitution model from IQ_TREE and parse it to RAxML-NG format


MODEL="$1"


MODEL="$(sed 's/K2P/K80/' <<<"$MODEL")"

MODEL="$(sed 's/TN/TN93/' <<<"$MODEL")"
MODEL="$(sed 's/TN93e/TN93ef/' <<<"$MODEL")" #  povodne v model bolo TNe, len predosly sed to zmenil na TN93e ;  TNe  je TN93ef
MODEL="$(sed 's/K3P/K81/' <<<"$MODEL")"
MODEL="$(sed 's/K81u/K81uf/' <<<"$MODEL")" #  to iste co vyssie
MODEL="$(sed 's/TPM2u/TPM2uf/' <<<"$MODEL")"
MODEL="$(sed 's/TPM3u/TPM3uf/' <<<"$MODEL")"

MODEL="$(sed 's/TIM+/TIM1uf+/' <<<"$MODEL")"
MODEL="$(sed 's/TIMe/TIM1/' <<<"$MODEL")" # opat to co vyssie, len naopak

MODEL="$(sed 's/TIM2/TIM2uf/' <<<"$MODEL")"
MODEL="$(sed 's/TIM2ufe/TIM2/' <<<"$MODEL")" # opat to co vyssie, len naopak

MODEL="$(sed 's/TIM3/TIM3uf/' <<<"$MODEL")"
MODEL="$(sed 's/TIM3ufe/TIM3/' <<<"$MODEL")" # opat to co vyssie, len naopak

MODEL="$(sed 's/TVMe/TVMef/' <<<"$MODEL")"  

echo "$MODEL" # odosli

exit



