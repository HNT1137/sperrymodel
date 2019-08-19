#Sperry Model, C++ Version
Ported from the original VBA implementation
Port and Readme author: Henry Todd (henry.todd@utah.edu)
(6/26/19) 0.23

------------

Introduction:

The model uses a stomatal gain vs. risk optimization (Sperry et al. 2017 -- see references below)
combined with a soil water budget to model plant responses to environmental conditions. Tracking
both the soil water content and the conductance loss due to cavitation allows the cumulative
effects of exposure to growing season weather to be modeled, including the effects of previous
drought events in the season.

Inputs are plant and site traits, considered to be static, and hourly weather drivers. The
outputs include net carbon assimilation (Anet), internal [CO2] (Ci), transpiration (E),
total evapotranspiration (ET), and element conductances (k) on an hourly and summary basis
(see details below).

------------

Basic usage instructions:

Provide plant and site traits through the parameters file, then supply hourly weather data to drive
model. The model will expect these files to be located in the current working directory. See the
included examples and the details below for formatting and units. These example files contain
all inputs necessary to test-run the model immediately after building.

To run, build and execute the model program (no command line arguments -- runs from files in
the working directory).
Building requires c++11, GNU example:
	> g++ -std=c++11
The -O3 and -ffast-math optimizations are recommended with GNU compilers:
	> g++ -std=c++11 -O3 -ffast-math
This version has also been tested with Visual Studio 2017's compiler with similar optimizations
(floating point mode fast, maximum optimization preferring speed).

The following files (included in this .zip) should be located in the working directory (normally
the same directory as the executable) before running:
	parameters.csv (plant, site params and program options)
	nametable.csv (maps parameter names to row/col locations in the parameters.csv sheet)
	dataset.csv (hourly weather drivers)
	dataheader.csv (a header row for the hourly data output)
	sumheader.csv (a header row for the summary data output)
	seasonlimits.csv (growing season limits, only required if using "sequential year mode"
	described below)

Upon completion, two output files are produced:
	-The "timesteps" output contains all of the hourly model outputs corresponding to the
	weather inputs.
	-The "summary" output includes various total values for each year (for example, net
	growing season productivity Anet, net transpiration E, etc.)

------------

Plant and Stand Parameters:
"parameters.csv"

Configure plant traits and other parameters in "parameters.csv" (expected input units are
indicated). There is also an Excel .xslx version of this file included which highlights the inputs
used in yellow and includes additional comments. The "parameters" sheet from this workbook can
be exported as "parameters.csv" after editing, or the parameters.csv file can be edited directly.

Noteworthy plant traits include (but are not limited to):
	-Whole plant kMax (saturated whole-plant conductance),
	-Percent of resistance in leaves (determines how tree conductance is partitioned to
	woody vs. leaf elements),
	-Vulnerability curves (in the form of Weibull curve B and C parameters),
	-Basal area/ground area (BA:GA, tree density or BAI),
	-Leaf area/basal area (together LA:BA and BA:GA determine LAI),
	-Leaf width,
	-Root Beta (controls rooting depth),
	-Maximum carboxylation rate at 25C, Vcmax25 (and associated maximum electron transport
	rate Jmax25, assumed to be Vmax25 * 1.67),

Important environment or site traits (again not comprehensive):
	-Ambient [CO2] Ca (input as ppm),
	-Soil hydraulic parameters,
	-Elevation,
	-Lat/lon,
	-Solar noon correction (offset between hour 12 in weather data and actual solar noon
	at this location)
	-Atmospheric "clear sky" transmittance (tau). Calibrates the amount of observed solar
	radiation considered to be "clear sky" (no clouds), generally between 0.6-0.75. See the
	equations in the "solarcalc" function if you would like to back-calculate transmittance
	from a observed clear sky data point.

-"parameters - inputs worksheet.xlsx" includes a "root and xylem worksheet" which can be used
to calculate Weibull B and C values for the vulnerability curve from P50 and P98 values. If
VC measurements are unavailable for certain elements of the plant other element curves can be
substituted. This sheet also includes calculations for converting root "beta" values to total
rooting depth in cm.

-The soil parameters we used for many soil types can be found in the "Common Soil Types" sheet
of "parameters - inputs worksheet.xlsx". Note: Currently only supports using the same soil type
for all active layers.

-Soil layers count (Default: 5) can be up to 5. A higher number of soil layers provides a more
robust soil water budget simulation, while fewer soil layers may improve performance slightly.

-Enabling ground water (Default: n) provides an unlimited source of water at a set potential and
distance below the root layers. This water will flow up into the soil layers, and potentially
allow layers to fill above field capacity (from the bottom layer up). When disabled (default),
the only sources of water input will be the initial fraction of field capacity and observed
rainfall (and any water over field capacity will become "drainage").

-Rain (Default: y) weather data rainfall will be ignored if disabled.

-Refilling (Default: n) allows trees to restore lost conductance, however the refilling model
is not sufficient to simulate authentic xylem refilling behavior and has not been thoroughly
tested in the current version of the code.

-Soil redistribution (Default: y) allows water to flow between soil layers

-Soil evaporation (Default: y) enables simulation of water evaporation from the surface soil layer.

-Use GS Data (Default: n) If enabled, multiple years will be run "sequentially" with on and
off seasons defined in seasonlimits.csv. See "Sequential year mode" for details. When disabled
(default), all weather timesteps provided are treated as part of the growing season and the
user is expected to truncate individual years to their start/end days. Water budget is reset
between years when disabled, treating years as totally independent.

-Autosave is always enabled regardless of the setting, as this version of the model has no
alternative output method. Output files will be generated in the working directory when the
run completes.

------------

Weather Data:
"dataset.csv"

-See example data for formatting. Weather drivers should be in hourly timesteps and can include
multiple years of data. Note that while year values are arbitrary, they should be sequential. For
example, if running data for the years 1997 and 2005 these should be numbered sequentially as
1 and 2 (or 1997 and 1998, etc).

-Inputs:
	-Year,
	-Julian Day (1-366),
	-Hour (0-23),
	-Obs. Solar (W m-2),
	-Rain (mm),
	-Wind (m s-1),
	-Tair (C),
	-Tsoil (C, if not available substitute air temp),
	-D (kPa)

------------

Outputs:

-Hourly Outputs (see dataheader.csv for full list):
	-Pressures (predawn soil layer pressures, sun and shade "mid-day" canopy pressures, MPa),
	-Water flows (mmol m-2s-1),
	-PS assimilation (A, umol s-1m-2 (leaf area)),
	-Gain-risk optimized stomatal conductance to water (Gw, mmol m-2s-1),
	-Element and whole plant conductances, (k, kghr-1m-2),
	-Water content and deltas (mm),
	-Ci

-Summary Outputs (per year, see sumheader.csv for full list):
	-Total Anet (mmol yr-1 m-2(leaf area)),
	-Total E (mm = mm3/mm2(ground area)),
	-Minimum whole plant conductance during the growing season (kghr-1m-2),
	-Percent Loss Conductance (PLC, percent, relative to a reference conductance at field
	capacity),
	-Mean Ci/Ca (+ A weighted Ci/Ca),
	-Water summary (start/end content, total growing season input (mm)).

------------

Sequential year processing: When running multiple years of data in a single dataset, the years
can be treated as entirely independent (the default) or can work from a continuous water budget.

*Independent year mode (default)
	-Set "Use GS Data" to "n" under "Program Options"
	-Use growing season trimmed data (see the example: "dataset.csv").
	-Ensure that the growing season limits are defined in "seasonlimits.csv"

The default setting is to reset the tree hydraulics and reset the soil water content to the
specified percent of field capacity every year. The years are completely independent, only run
in a single dataset for convenience.

In this mode, weather data should be trimmed to only the growing season days as in the included
dataset.csv (so that the last day of one growing season is followed immediately by the first day
of the next). Year values are used for output, and each must be unique and optimally sequential
(to determine when new years begin), but the values are otherwise unused by the model and
thus do not need to be meaningful values. When running in this mode the growing season limits
(seasonlimits.csv) will not be used; All days in the dataset will be considered to be in the
growing season.

-------------

*Sequential year mode:
	-Set "Use GS Data" to "y" under "Program Options"
	-Use full-year data (see the example: "dataset - full year example.csv")
	-Ensure that the growing season limits are defined in "seasonlimits.csv"
		--Note that the year values in "seasonlimits.csv" are for reference only;
		The first row of start/end days will be used for the first year of data, etc.

In this mode, plant hydraulics will reset between seasons and plant transpiration/productivity will
be disabled during the off-season, but soil water budget will continue to be computed. Soil may
or may not be refilled to field capacity depending on the availability of off-season precipitation.

Note that soil surface evaporation will also be disabled during the off-season. This is not
particularly realistic, but the functionality was intended to answer the question: Is there at
minimum enough recorded rainfall to refill the soil? A more robust off-season water simulation
would require additional data (snow pack) and simulation of soil behavior under snow and is
not provided here.

-------------

References:

Describing the gain/risk algorithm used in the model:

	-Sperry JS, Venturas MD, Anderegg WRL, Mencucinni M, Mackay DS, Wang Y, Love DM (2017)
	Predicting stomatal responses to the environment from the optimization of photosynthetic
	gain and hydraulic cost. Plant Cell and Environment 40: 816-830

Describing the original hydraulic model the gain-risk optimization was based on:

	-Sperry JS, Love DM (2015) Tansley Review: What plant hydraulics can tell us about
	plant responses to climate-change droughts. New Phytologist 207: 14-17

	-Sperry JS, Wang Y, Wolfe BT, Mackay DS, Anderegg WRL, McDowell NG, Pockman WT (2016)
	Pragmatic hydraulic theory predicts stomatal responses to climatic water deficits. New
	Phytologist 212: 577-589

(Full-text PDFs available at http://sperry.biology.utah.edu/publications/)

-------------

Contact:

For specific questions about this C++ version of the model, contact:
Henry Todd
henry.todd@utah.edu
hnt1137@gmail.com
