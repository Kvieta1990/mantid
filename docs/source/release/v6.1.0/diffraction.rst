===================
Diffraction Changes
===================

.. contents:: Table of Contents
   :local:

.. warning:: **Developers:** Sort changes under appropriate heading
    putting new features at the top of the section, followed by
    improvements, followed by bug fixes.

Powder Diffraction
------------------
- New algorithm :ref:`RebinRagged <algm-RebinRagged>` which can rebin a workspace with different binning parameters for each spectrum
- New options, including three "cache directory" and one "clean cache" in the Advanced Setup tab of the SNS Powder Reduction interface
- New caching feature is added to :ref:`SNSPowderReduction <algm-SNSPowderReduction>` to speed up calculation using same sample and container.
- Use caching for Va in :ref:`SNSPowderReduction <algm-SNSPowderReduction>`.
- Improve algorithm :ref:`FitPeaks <algm-FitPeaks>` to enable it to fit with multiple peaks in same spectrum with Back-to-back Exponential function starting from user specified parameters.
- :ref:`SNSPowderReduction <algm-SNSPowderReduction>` has additional property, ``DeltaRagged``, which allows using :ref:`RebinRagged <algm-RebinRagged>` to bin each spectrum differently.

Engineering Diffraction
-----------------------

Single Crystal Diffraction
--------------------------
- New version of algorithm :ref:`SCDCalibratePanels <algm-SCDCalibratePanels-v2>` provides more accurate calibration results for CORELLI instrument.
- Modified some logs in output workspace from :ref:`LoadWANDSCD <algm-LoadWANDSCD>` to be TimeSeriesProperty so they work with :ref:`SetGoniometer <algm-SetGoniometer>`.

Improvements
############
- :ref:`IntegratePeaksMD <algm-IntegratePeaksMD>` now allows ellipsoidal shapes to be manually defined for the PeakRadius and Background radii options.
- :ref:`SNSPowderReduction <algm-SNSPowderReduction>` now check if previous container is created using the same method before reusing it.

:ref:`Release 6.1.0 <v6.1.0>`
