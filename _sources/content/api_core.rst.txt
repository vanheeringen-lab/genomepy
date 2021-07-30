Python API documentation (core)
===============================

.. currentmodule:: genomepy

This page described the core genomepy functionality.
These classes and functions can be found on the top level of the genomepy module (e.g. :code:`genomepy.search`),
and are made available when running :code:`from genomepy import *` (we won't judge you).

Additional functions that do not fit the core functionality, but we feel are still pretty cool, are also described.

Finding genomic data
--------------------

When looking to download a new genome/gene annotation, your first step would be :code:`genomepy.search`.
This function will check either one, or all, providers.
Advanced users may want to specify a provider for their search to speed up the process.
To see which providers are available, use :code:`genomepy.list_providers`,
or :code:`genomepy.list_online_providers`:

------------

.. autofunction:: list_providers
   :noindex:

.. autofunction:: list_online_providers
   :noindex:

.. autofunction:: search
   :noindex:

------------

If you have no idea what you are looking for, you could even check out *all* available genomes.
Be warned, :code:`genomepy.list_available_genomes` is like watching the Star Wars title crawl.

------------

.. autofunction:: list_available_genomes
   :noindex:

------------

If we search for homo sapiens for instance, we find that :code:`GRCh3.p13` and :code:`hg38` are the latest versions.
These names describe the same genome, but different :code:`assemblies`, with differences between them.

One of these differences is the quality of the gene annotation.
Next, we can inspect these with :code:`genomepy.head_annotations`:

------------

.. autofunction:: head_annotations
   :noindex:

------------

Installing genomic data
------------------------

Now that you have seen whats available, its time to download a genome.
The default parameter for :code:`genomepy.install_genome` are optimized for sequence alignment and gene counting,
but you have full control over them, so have a look!

genomepy won't overwrite any files you already downloaded (unless specified),
but you can review your local genomes with :code:`genomepy.list_installed_genomes`.

------------

.. autofunction:: install_genome
   :noindex:

.. autofunction:: list_installed_genomes
   :noindex:

------------

If you want to download a sequence blacklist, or create an aligner index, you might wanna look at plugins!
Don't worry, you can rerun the :code:`genome.install_genome` command, and genomepy will only run the new parts.

------------

.. autofunction:: manage_plugins
   :noindex:

------------

The genome and gene annotations were installed in the genomes directory (unless specified otherwise).
If you have a specific location in mind, you could set this as default in the genomepy config.
To find and inspect it, use :code:`genomepy.manage_config`:

------------

.. autofunction:: manage_config
   :noindex:

------------

Errors
------

Did something go wrong? Oh noes!
If the problem persists, clear the genomepy cache with :code:`genomepy.clean`, and try again.

------------

.. autofunction:: clean
   :noindex:

------------

Using a genome
--------------

Alright, you've got the goods!
You can browse the genome's sequences and metadata with the :code:`genomepy.Genome` class.
This class builds on the :code:`pyfaidx.Fasta` class to also provide you with several options
to get specific sequences from your genome, and save these to file.

------------

.. autoclass:: Genome
   :noindex:
   :members:
   :inherited-members:

   .. rubric:: Methods
   .. autosummary::
      ~Genome.close
      ~Genome.get_random_sequences
      ~Genome.get_seq
      ~Genome.get_spliced_seq
      ~Genome.items
      ~Genome.keys
      ~Genome.track2fasta
      ~Genome.values

   .. rubric:: Attributes
   .. autosummary::
      ~Genome.gaps
      ~Genome.plugin
      ~Genome.sizes
      ~Genome.genomes_dir
      ~Genome.name
      ~Genome.genome_file
      ~Genome.genome_dir
      ~Genome.index_file
      ~Genome.sizes_file
      ~Genome.gaps_file
      ~Genome.annotation_gtf_file
      ~Genome.annotation_bed_file
      ~Genome.readme_file

..
   Including the autosummary file works, but there's a tonne of warnings.
   Plus, the class (& contents) are all marked as genomepy.genome.Genome...
   The code above is mostly copy-pasted from this file now.
   .. include:: ../_autosummary/genomepy.genome.Genome.rst

------------

You can obtain genomic sequences from a wide variety of inputs with :code:`as_seqdict`.
To use the function, it must be explicitly imported with :code:`from genomepy.seq import as_seqdict`.

------------

.. autofunction:: genomepy.seq.as_seqdict
   :noindex:

------------

A non-core function worth mentioning is :code:`genomepy.files.filter_fasta`,
for if you wish to filter a fasta file by chromosome name using regex,
but want the output straight to (another) fasta file.

------------

.. autofunction:: genomepy.files.filter_fasta
   :noindex:

------------

Using a gene annotation
-----------------------

Similarly, the :code:`genomepy.Annotation` class helps you get the genes in check.
This class returns a number of neat pandas dataframes, such as the :code:`named_gtf`,
or an annotation with the gene or chromosome names remapped to another type.
Remapping gene names to another type is also with the Annotation class with :code:`Annotation.map_genes`,
but also separately with :code:`genomepy.query_mygene`, as it's just so damn useful.

------------

.. autoclass:: Annotation
   :members:
   :noindex:

.. autofunction:: query_mygene
   :noindex:

------------

Another non-core function worth mentioning is :code:`genomepy.annotation.filter_regex`,
which allows you to filter a dataframe by any columns using regex.

------------

.. autofunction:: genomepy.annotation.filter_regex
   :noindex:
