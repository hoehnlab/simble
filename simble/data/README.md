# Data Files

This repository contains data files that were obtained from various sources. Below is a list of the sources for each file:

- `hh_sf5.csv` was generated in R with [shazam](https://shazam.readthedocs.io/en/stable/) using
`shazam::HH_S5F@mutability` [HH_S5F][HH_S5F][^HH_S5F]. This file represents the normalized rates of a given 5-mer being mutated in heavy chain sequences from Homo sapiens.

- `hh_sf5_substitution.csv` was generated in R with [shazam](https://shazam.readthedocs.io/en/stable/) using
`shazam::HH_S5F@substitution` [HH_S5F][HH_S5F][^HH_S5F]. This file represents the normalized rates of the center nucleotide of a given 5-mer mutating to a different nucleotide in heavy chain sequences from Homo sapiens.

- `hkl_sf5.csv`: was generated in R with [shazam](https://shazam.readthedocs.io/en/stable/) using
`shazam::HKL_S5F@mutability` [HKL_S5F][HKL_S5F][^HKL_S5F]. This file represents the normalized rates of a given 5-mer being mutated in kappa and lambda light chain sequences from Homo sapiens.

- `hkl_sf5_substitution.csv`: was generated in R with [shazam](https://shazam.readthedocs.io/en/stable/) using
`shazam::HKL_S5F@msubstitution` [HKL_S5F][HKL_S5F][^HKL_S5F]. his file represents the normalized rates of the center nucleotide of a given 5-mer mutating to a different nucleotide in kappa and lambda light chain sequences from Homo sapiens.

- `naive_heavy_chain.tsv`: correspondence from Kenneth B. Hoehn

- `naive_light_chain.tsv`: correspondence from Kenneth B. Hoehn

- `naive_pairs_filtered.csv`: The contents of this file were generated from `naive_heavy_chain.tsv` and `naive_light_chain.tsv` using the custom script `clean_naive.py` in that also appears in this folder.


Please note that the data files are provided for informational purposes only and should be used responsibly and in accordance with the terms and conditions of the respective sources.


<!-- If you have any questions or concerns regarding the data files, please contact us at [contact@example.com](mailto:contact@example.com). -->


[HH_S5F]: https://doi.org/10.3389/fimmu.2013.00358 "Yaari G, Vander Heiden J, Uduman M, Gadala-Maria D, Gupta N, Stern J, O’Connor K, Hafler D, Lasserson U, Vigneault F, Kleinstein S (2013). “Models of somatic hypermutation targeting and substitution based on synonymous mutations from high-throughput immunoglobulin sequencing data.” Frontiers in Immunology, 4(358), 1-11. doi:10.3389/fimmu.2013.00358 https://doi.org/10.3389/fimmu.2013.00358."

[HKL_S5F]: https://doi.org/10.4049/jimmunol.1502263 "Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O’Connor K, Vigneault F, Shlomchik M, Kleinstein S (2016). “A Model of Somatic Hypermutation Targeting in Mice Based on High-Throughput Ig Sequencing Data.” The Journal of Immunology, 197(9), 3566-3574. doi:10.4049/jimmunol.1502263 https://doi.org/10.4049/jimmunol.1502263."


[^HH_S5F]: Yaari G, Vander Heiden J, Uduman M, Gadala-Maria D, Gupta N, Stern J, O’Connor K, Hafler D, Lasserson U, Vigneault F, Kleinstein S (2013). “Models of somatic hypermutation targeting and substitution based on synonymous mutations from high-throughput immunoglobulin sequencing data.” Frontiers in Immunology, 4(358), 1-11. doi:10.3389/fimmu.2013.00358 [https://doi.org/10.3389/fimmu.2013.00358](https://doi.org/10.3389/fimmu.2013.00358).

[^HKL_S5F]: Cui A, Di Niro R, Vander Heiden J, Briggs A, Adams K, Gilbert T, O’Connor K, Vigneault F, Shlomchik M, Kleinstein S (2016). “A Model of Somatic Hypermutation Targeting in Mice Based on High-Throughput Ig Sequencing Data.” The Journal of Immunology, 197(9), 3566-3574. doi:10.4049/jimmunol.1502263 [https://doi.org/10.4049/jimmunol.1502263](https://doi.org/10.4049/jimmunol.1502263).