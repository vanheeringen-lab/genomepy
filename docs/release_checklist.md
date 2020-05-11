# Release checklist

1. Create release candidate with `git flow`:

```
$ git flow release start ${new_version} 
```

2. Update version in `genomepy/__about__.py`

3. Make sure `CHANGELOG.md` is up-to-date.

* you can commit, but do not need to push to remote


4. Make sure all tests pass.

5. Check if release works on pypi:

```
python setup.py sdist bdist_wheel
twine upload --repository-url https://test.pypi.org/legacy/ dist/genomepy-${version}*

pip install --extra-index-url https://test.pypi.org/simple/ genomepy==${version}
genomepy search xenopus_tropicalis
```

6. Finish the release:

```
git flow release finish ${new_version}
```

7. Push everything to github, including tags:

```
git push --follow-tags origin develop
```

8. Pull into master
  
9. Upload to pypi:

```
python setup.py sdist bdist_wheel
twine upload dist/genomepy-${version}*
```

10. Create release on github (if it not already exists)

* Update release with CHANGELOG information from the latest version
* Download the tarball from the github release (`.tar.gz`). 
* Attach downloaded tarball to release as binary (this way the download count get tracked).


11. Update bioconda package

* fork bioconda/bioconda-recipes
* follow the steps in the [docs](https://bioconda.github.io/contributor/workflow.html)
* get the hash from the downloaded tarbal using `sha256sum *genomepy-${version}.tar.gz`
* update the [yaml file](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/genomepy/meta.yaml) locally. 
* push to a new branch on the fork
* start a PR
