.. NEMO documentation master file, created by
   sphinx-quickstart on Tue Oct  9 18:14:01 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#################
Quick Start Guide
#################

..
	A hidden .global.rst should be included in every subfiles with `include` directive
   It contains a list of common URL links 
     
.. include:: .global.rst

.. include:: readme.rst

Summary
=======

.. toctree::
	:maxdepth: 1
	:titlesonly:

	release_notes.rst
	install.rst
	reference_configurations.rst
	test_cases.rst
	setup_configuration.rst
	interfacing_options.rst
	definitions.rst

..
   For headings markup, this convention is recommended from Python’s Style Guide
	# with overline, for parts
	* with overline, for chapters
	=, for sections
	-, for subsections
	^, for subsubsections
	", for paragraphs

..
	Indices and tables
	==================
	* :ref:`genindex`
	* :ref:`modindex`
	* :ref:`search`
