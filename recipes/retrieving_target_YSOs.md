[retrieving_target_YSOs.md](https://github.com/user-attachments/files/22308720/retrieving_target_YSOs.md)
\# `retrieving_target_YSOs`


\*\*Authors: William Salazar, Eliza Fabian, \& Austin Jones\*\*


\## Description

This recipe guides you through retrieving YSOs from the SPICY catalog ([Kuhn et al. 2021](https://iopscience.iop.org/article/10.3847/1538-4365/abe465)). We seek to first retrieve the infrared magnitudes provides by the SPICY catalog. Filters we use include Spitzer's [IRAC bands](https://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=Spitzer&asttype=) and 2MASS [J,H,K bands](https://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=2MASS&asttype=). We will retrieve files needed for both the infrared only and with [Gaia bands](https://svo2.cab.inta-csic.es/theory/fps/index.php?mode=browse&gname=GAIA&gname2=GAIA3&asttype=) runs of the pipeline. At the time of writing this, the latest Gaia release is DR3.


## INITIAL SETUP

For this procedure, we solely need TOPCAT to slice the fits file from SPICY. If on a Windows operating system, you will also need Java. For TOPCAT, you can install it at https://www.star.bris.ac.uk/~mbt/topcat/ 

You will also need the SPICY catalog dataset, which will be provided to you or you can retrieve it from https://irsa.ipac.caltech.edu/data/SPITZER/GLIMPSE/overview.html

OPTIONAL: If you would like to visualize the data, you can use DS9, which can be downloaded at https://sites.google.com/cfa.harvard.edu/saoimageds9

## Process

# <font color='red'>IR-Only Dataset</font>

In TOPCAT, open the fits file to access the SPICY catalog. We also need galactic latitude (b) and longitude (l) coordinates of your region for the cutout.

To create a cutout, we go to "Display row subsets" (it looks like blue oval with a red oval inside) and click on the Plus sign to add in a new subset. It must be in the format 'l<l_max && l>l_min && b<b_max && b>b_min'. An example of a cutout would be:

- l<13.15 && l>12.4 && b<0.35 && b>-0.6

Then go to "Display column metadata" (icon with a blue row separated from the grey rows). We only keep the following table headers: 

- l, b, Spitzer, mag3_6, e_mag3_6, mag4_5, e_mag4_5, mag5_8, e_mag5_8, mag8_0, e_mag8_0, j_synth,h_synth,k_synth 

To save the subset, make sure under "Row Subsets" in the main console you have the subset name in place. Then on the Floppy Disk, save table. Save as the name of your region etc. name_of_region.csv . Now you have your file to run through the pipeline below.

Then, if you wish to save the session, click on the Floppy Disk and go under "Save Session". Then name the session and save.



# <font color='blue'>With Gaia Dataset</font>

To get the Gaia dataset, first we must go to [Gaia Archive](https://gea.esac.esa.int/archive/), Then, find the SEARCH tab -> Advanced (ADQL). In the box type the following:

SELECT
   TOP 1000000
   *
   FROM gaiadr3.gaia_source
   WHERE
	l BETWEEN l_min AND l_max
  	AND
  	b BETWEEN b_min AND b_max

Next, Download format: CSV (dropdown at bottom of the page). Download by clicking on the button that looks like a stack of pancakes with a + sign. You can replace ra and dec with l and b and your galactic longitude and latitude numbers.


Using the csv file you downloaded from the Gaia Archive, open it in TOPCAT. We want to first match and combine the table from the IR only dataset and Gaia. Using the "Create new table by matching rows in two existing table" button (matchsticks icon), we match on the GAIADR2 column in the IR data and source_id in the Gaia dataset. Use the "Exact Value" option on the "Algorithm" dropbox (can also do this using the "SKY" option by matching the l and b columns and using 1 arcsecond or less as the matching criteria.) For this table, you want to save the following table headers: 

- Spitzer, phot_bp_mean_flux_over_error, phot_bp_mean_mag, phot_rp_mean_flux_over_error, phot_rp_mean_mag

IMPORTANT: Save the "over_error" columns (as titled above), not the "error" columns, as those are not the true errors.

NOTE: The Gaia G band is not saved since it overlaps with the other two bands.


## NEXT STEPS:

Now proceed to the first pipeline recipe 'data_format_pipeline'.
