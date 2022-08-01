:orphan:

# Release checklist

1. Make sure all tests pass.

    ```shell
    mamba env update -f environment.yml
    pytest -vvv
    ```

2. Create release candidate with `git flow`:

    ```shell
    new_version=0.0.0
    echo ${new_version}
    
    git flow release start ${new_version}
    ```

3. Update version in `genomepy/__about__.py`

4. Update the `CHANGELOG.md`

    * add the new version & date to the header
    * link to the diff in the footer
    * add & commit the changes, but do not push

5. Check if release works on pypi:

    ```shell
    python -m build

    # twine must be up to date (3.3.0 works). System installed twine can interfere.
    twine upload --repository-url https://test.pypi.org/legacy/ dist/genomepy-${new_version}*

    pip uninstall genomepy

    # the \ is to escape the ==, so the variable ${new_version} can be called
    pip install --extra-index-url https://test.pypi.org/simple/ genomepy\==${new_version}

    # tests
    genomepy --version
    genomepy --help
    genomepy install --help
    genomepy clean
    genomepy search xenopus_tropicalis
    genomepy annotation hg38
    genomepy annotation GRCh38.p13
    genomepy install -af -p gencode GRCm39
    genomepy install -af -p ensembl TAIR10
    genomepy install -af -p ucsc sacCer3 --UCSC-annotation ensGene
    genomepy install -af -p ncbi ASM2732v1
    genomepy install -af -p url -l url_test https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz --URL-to-annotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz
    genomepy install -af -p local -l local_test ~/.local/share/genomes/TAIR10/TAIR10.fa --Local-path-to-annotation ~/.local/share/genomes/TAIR10/TAIR10.annotation.gtf
    
    ```

6. Finish the release:

    ```shell
    git flow release finish ${new_version}
    ```

7. Push everything to github, including tags:

    ```shell
    git push --follow-tags origin develop master
    ```

8. Upload to pypi:

    ```shell
    python -m build
    twine upload dist/genomepy-${new_version}*
    ```

9. Create release on github (if it not already exists)

    * Update release with CHANGELOG information from the latest version
    * Download the tarball from the github release (`.tar.gz`).
    * Attach downloaded tarball to release as binary (this way the download count get tracked).

10a. Update bioconda package

    * wait for the bioconda bot to create a PR
    * update dependencies in the bioconda recipe.yaml if needed
    * approve the PR
    * comment: @bioconda-bot please merge

10b. Update bioconda package

    * fork bioconda/bioconda-recipes
    * follow the steps in the [docs](https://bioconda.github.io/contributor/workflow.html)
    * get the hash from the downloaded tarbal using `sha256sum *genomepy-${new_version}.tar.gz`
    * update the [yaml file](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/genomepy/meta.yaml) locally.
    * push to a new branch on the fork
    * start a PR
