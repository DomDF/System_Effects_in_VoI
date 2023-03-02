#
# "System Effects in Identifying Risk-Optimal Data Requirements for Digital Twins of Structures"
#  Paper submitted to the journal of Reliability Engineering & System Safety
#  Domenic Di Francesco, PhD, CEng (MIMechE)
#  The Alan Turing Institute, University of Cambridge
#  
#  Accompanying Julia code to set-up and run calculations
#

######################################################
#
# Loading libraries
#
######################################################

# For describing probabilistic models
using Distributions, Turing, Random, LatinHypercubeSampling, Copulas
# For describing and solving decision problem
using JuMP, HiGHS, DecisionProgramming, LinearAlgebra
# For working with data
using CSV, DataFrames, DataFramesMeta