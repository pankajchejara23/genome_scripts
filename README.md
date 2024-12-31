# **Genome Data Processing Scripts**

This repository contains a set of scripts designed to automate tasks in genome sequence analysis. The scripts help streamline the workflow for processing raw genomic data, enabling efficient downstream analysis.

---

## **Table of Contents**
1. [Overview](#overview)
2. [Scripts & Their Functionality](#scripts--their-functionality)
3. [Prerequisites](#prerequisites)
4. [Getting Started](#getting-started)
5. [License](#license)

---

## **Overview**

The scripts in this repository automate several key steps in genome sequence analysis, including downloading raw sequencing data, applying quality control measures, and preparing data for downstream analysis. These tools are especially useful for bioinformatics researchers working with high-throughput sequencing data.

---

## **Scripts & Their Functionality**

Below is a list of the available scripts along with their intended functionality:

| **Script**         | **Functionality**                                                                | **Prerequisite** |
|--------------------|---------------------------------------------------------------------------------|------------------|
| `ncbi_downloader`   | Downloads raw read sequences and metadata from the NCBI SRA for a specified project. | **SRA Toolkit**  |
| `apply_trimmomatic` | Applies Trimmomatic for trimming paired-end read sequences.                      | **Trimmomatic**  |
| `create_manifest`   | Creates a manifest file for QIIME2 import, facilitating analysis in QIIME2.       | --               |

---

## **Prerequisites**

To use these scripts, you will need to have the following tools installed:

1. **SRA Toolkit** – Required for downloading sequence data from NCBI SRA.
   - Installation: [SRA Toolkit installation guide](https://github.com/ncbi/sra-tools/wiki/Downloads)
2. **Trimmomatic** – Required for quality trimming of read sequences.
   - Installation: [Trimmomatic download page](http://www.usadellab.org/cms/?page=trimmomatic)
3. **QIIME2** – Required for downstream analysis in QIIME2 (for the `create_manifest` script).
   - Installation: [QIIME2 installation guide](https://docs.qiime2.org/2023.2/install/)

---

## **Getting Started**

To get started with these scripts, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/genome-data-processing-scripts.git
   cd genome-data-processing-scripts

2. Make sure you have all the prerequisites installed as listed above.

3. Run individual scripts as needed

## **License**

This repository is licensed under the MIT License.