Coupling with other models (OASIS, SAS, ...)
============================================

.. include:: .global.rst

NEMO currently exploits OASIS-3-MCT to implement a generalised coupled interface
(`Coupled Formulation <http://forge.ipsl.jussieu.fr/nemo/doxygen/node50.html?doc=NEMO>`_).
It can be used to interface with most of the European atmospheric GCM (ARPEGE, ECHAM, ECMWF, Ha- dAM, HadGAM, LMDz),
as well as to WRF (Weather Research and Forecasting Model), and to implement the coupling of
two independent NEMO components, ocean on one hand and sea-ice plus other surface processes on the other hand
(`Standalone Surface Module - SAS <http://forge.ipsl.jussieu.fr/nemo/doxygen/node46.html?doc=NEMO>`_).

To enable the OASIS interface the required compilation key is ``key_oasis3``.
The parameters to set are in sections ``namsbc_cpl`` and in case of using of SAS also in section ``namsbc_sas``.
