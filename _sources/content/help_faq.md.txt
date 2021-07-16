# Frequently Asked Questions

* Q: genomepy sais a provider is offline!

  A: In rare occasions providers can be down for several hours.
  However, in most situations the provider was just busy, and will be available again in a minute or so.

* Q: There is something wrong with a provider!

  A: If the provider was offline when genomepy downloaded their database the local data may be incorrect.
  To fix this, simply run `genomepy clean` on the command line, or run `genomepy.clean()` in Python.

* Q: which genome/gene annotation should I use?

  A: each provider has its pros and cons:
    * Ensembl has excellent gene annotations, but their chromosome names can cause issues with some tools.
    * UCSC has an excellent genome browser, but their gene annotations vary in format.
    * NCBI allows public submissions, and so has the latest versions, although not always complete or error free.

  Use `genomepy search` to see your options, and `genomepy annotation` to check the quality of the gene annotation(s).

* Q: genomepy search showed me a genome, but the download didn't work!

  A: It is possible the provider made a typo in the genome/annotation URL.
  If this is *not* the case, please make an issue on [our github page](https://github.com/vanheeringen-lab/genomepy/issues)!

* Q: genomepy keeps using a plugin, how do I turn those off again?

  A: On the command line, you can use `genomepy plugin --help` to see the options.

* Q: my genomepy config was corrupted, help!

  A: you can create a new one with `genomepy config generate` on command line,
  or `genomepy.manage_config("generate")` in Python.
