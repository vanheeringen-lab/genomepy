# Release checklist

1. Create release candidate with `git flow`.

```
$ git flow release start ${new_version} 
```

2. Update version in `genomepy/__about__.py`

3. Make sure `CHANGELOG.md` is up-to-date.

4. Check if release works on pypi

```
python setup.py sdist bdist_wheel
twine upload --repository-url https://test.pypi.org/legacy dist/genomepy-${version}*
```

5. Finish the release

```
git flow release finish ${new_version}
```

6. Push everything to github, including tags


