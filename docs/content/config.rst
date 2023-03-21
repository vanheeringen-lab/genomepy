=============
Configuration
=============

Genomepy uses many configurable options.
The default setting for several options can be overwritten using a config file.

You can create new config file with :code:`genomepy config generate`,
or find the location of your current config file with :code:`genomepy config file`.

Configurable options
--------------------

The config file uses the YAML format, wherein each configurable option is given as a :code:`key: value` pair,
and anything after the comment symbol :code:`#` is ignored.

These keys are currently supported:

* :code:`bgzip`: determines if newly installed assembly data is compressed or not.

    * Options: :code:`true` or :code:`false`.

    * Default: :code:`false`

* :code:`cache_exp_genomes`: expiration time (in seconds) for the cache used in :code:`genomepy search`.
  (Re)building this cache is slow, but must be done periodically to get the latest assemblies.
  This setting does not affect your installed assembly data.

    * Options: an integer or scientific number, or :code:`None` for infinite.

    * Default: 6.048e5 (1 week)

* :code:`cache_exp_other`: expiration time (in seconds) for the cache used in short term activities (e.g. mygene.info queries).
  This setting does not affect your installed assembly data.

    * Options: an integer or scientific number, or :code:`None` for infinite.

    * Default: 3.6e3 (1 hour)

* :code:`genomes_dir`: the path where assembly data will be installed.

    * Default: :code:`~/.local/share/genomes/`

* :code:`plugin`: the list of currently active plugins.
  See command :code:`genomepy plugin show` for options.
  :code:`genomepy plugin` can also be used to easily change this list.
