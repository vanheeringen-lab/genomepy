:orphan:

# Release checklist

1. Make sure all tests pass.

    `pytest -vv --disable-pytest-warnings`

2. Create release candidate with `git flow`:

    ```
    new_version=0.0.0
    git flow release start ${new_version}
    ```

3. Update version in `genomepy/__about__.py`

4. Update the `CHANGELOG.md`

    * add the new version & date to the header
    * link to the diff in the footer
    * add & commit the changes, but do not push

5. Check if release works on pypi:

    ```shell
    python setup.py sdist bdist_wheel

    # twine must be up to date (3.3.0 works). System installed twine can interfere.
    twine upload --repository-url https://test.pypi.org/legacy/ dist/genomepy-${new_version}*

    python setup.py develop --uninstall
   
    # the \ is to escape the ==, so the variable ${new_version} can be called
    pip install --extra-index-url https://test.pypi.org/simple/ genomepy\==${new_version}

    # tests
    genomepy --version
    genomepy --help
    genomepy install --help
    genomepy clean
    genomepy search xenopus_tropicalis
    genomepy install -af -p ensembl TAIR10
    genomepy install -af -p ucsc sacCer3
    genomepy install -af -p ncbi ASM2732v1
    genomepy install -af -p url  -l url_test  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.fna.gz --URL-to-annotation https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/027/325/GCF_000027325.1_ASM2732v1/GCF_000027325.1_ASM2732v1_genomic.gff.gz
    ```

6. Finish the release:

    `git flow release finish ${new_version}`

7. Push everything to github, including tags:

    `git push --follow-tags origin develop`

8. Pull into master
  
9. Upload to pypi:

    ```
    python setup.py sdist bdist_wheel
    twine upload dist/genomepy-${new_version}*
    ```

10. Create release on github (if it not already exists)

    * Update release with CHANGELOG information from the latest version
    * Download the tarball from the github release (`.tar.gz`).
    * Attach downloaded tarball to release as binary (this way the download count get tracked).

11a. Update bioconda package

    * wait for the bioconda bot to create a PR
    * update dependencies in the bioconda recipe.yaml if needed
    * approve the PR
    * comment: @bioconda-bot please merge

11b. Update bioconda package

    * fork bioconda/bioconda-recipes
    * follow the steps in the [docs](https://bioconda.github.io/contributor/workflow.html)
    * get the hash from the downloaded tarbal using `sha256sum *genomepy-${new_version}.tar.gz`
    * update the [yaml file](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/genomepy/meta.yaml) locally.
    * push to a new branch on the fork
    * start a PR
