#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
python find_enrichment.py study.file population.file gene-association.file
This program returns P-values for functional enrichment in a cluster of study
genes using Fisher's exact test, and corrected for multiple testing (including
Bonferroni, Holm, Sidak, and false discovery rate).
About significance cutoff:
--alpha: test-wise alpha; for each GO term, what significance level to apply
        (most often you don't need to change this other than 0.05 or 0.01)
--pval: experiment-wise alpha; for the entire experiment, what significance
        level to apply after Bonferroni correction
"""

__copyright__ = "Copyright (C) 2010-present, H Tang et al. All rights reserved."
__author__ = "various"

import sys
import os
from goatools.cli.find_enrichment import GoeaCliArgs
from goatools.cli.find_enrichment import GoeaCliFnc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def main():
    """Run gene enrichment analysis."""
    # Load study, population, associations, and GoDag. Run GOEA.
    obj = GoeaCliFnc(GoeaCliArgs().args)
    # Reduce results to significant results (pval<value)
    results_specified = obj.get_results()
    # Print results in a flat list
    obj.prt_results(results_specified)

if __name__ == "__main__":
    main()

# Copyright (C) 2010-present, H Tang et al. All rights reserved.
