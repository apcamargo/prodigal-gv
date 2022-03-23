# `prodigal-gv`

A fork of [Prodigal](https://github.com/hyattpd/Prodigal) meant to improve gene calling for giant viruses and viruses that use alternative genetic codes. It includes six models added to the metagenome mode:
* *Acanthamoeba polyphaga mimivirus*
* *Paramecium bursaria* Chlorella virus
* *Acanthocystis turfacea* Chlorella virus
* [VirSorter2](https://github.com/jiarong/VirSorter2)'s NCLDV gene model
* [Topaz (genetic code 15)](https://www.biorxiv.org/content/10.1101/2021.08.26.457843v1.full)
* [Agate (genetic code 15)](https://www.biorxiv.org/content/10.1101/2021.08.26.457843v1.full)

## Installation

`prodigal-gv` can be installed via [Conda/Mamba](https://anaconda.org/bioconda/prodigal-gv), pre-build binaries ([Linux](https://github.com/apcamargo/prodigal-gv/releases/download/2.7.0/prodigal-gv-linux) and [macOS](https://github.com/apcamargo/prodigal-gv/releases/download/2.7.0/prodigal-gv-macos)), or built from scratch:

```
# Conda/Mamba:
conda install -c bioconda prodigal-gv

# Pre-built binary (Linux):
curl -L https://github.com/apcamargo/prodigal-gv/releases/download/2.7.0/prodigal-gv-linux -o prodigal-gv
chmod +x prodigal-gv

# Pre-built binary (macOS):
curl -L https://github.com/apcamargo/prodigal-gv/releases/download/2.7.0/prodigal-gv-macos -o prodigal-gv
chmod +x prodigal-gv

# Buld from scratch
git clone https://github.com/apcamargo/prodigal-gv.git
cd prodigal-gv
make
```

## Usage

```bash
prodigal-gv -p meta -i genome.fna -a proteins.faa > /dev/null 2>&1
```

---

Prodigal was written by [Doug Hyatt](https://github.com/hyattpd/) and its usage should be acknowledged.

> Hyatt, D., Chen, G.-L., LoCascio, P.F., Land, M.L., Larimer, F.W., and Hauser, L.J. (2010). [Prodigal: prokaryotic gene recognition and translation initiation site identification](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-119). BMC Bioinformatics *11*, 119.
