Visualisation for Blast XML outputs under Galaxy
================================================


Our project
-----------

Blastvisu is a project support by [ABiMS](http://abims.sb-roscoff.fr/)
It aim to provide the equivalent of NCBI-W3BLAST and [NCBI-BLAST](http://blast.ncbi.nlm.nih.gov/Blast.cgi) on Galaxy.


Installation
------------

``` {.bash}
git clone https://github.com/lecorguille/galaxy-blastvisu.git
cp -r galaxy-blastvisu $GALAXY_ROOT/config/plugins/visualizations/blastvisu
```

Restart Galaxy



Blast
------

The Basic Local Alignment Search Tool (BLAST) finds regions of local similarity between sequences. The program compares nucleotide or protein sequences to sequence databases and calculates the statistical significance of matches. BLAST can be used to infer functional and evolutionary relationships between sequences as well as help identify members of gene families. 

Source: [http://blast.ncbi.nlm.nih.gov/Blast.cgi](http://blast.ncbi.nlm.nih.gov/Blast.cgi)


Galaxy
------
Galaxy is an open, web-based platform for data intensive biomedical research. Whether on the free public server or your own instance, you can perform, reproduce, and share complete analyses. 

Homepage: [https://galaxyproject.org/](https://galaxyproject.org/)


Galaxy_Blast
------------

The visualisation is design to display the XML format which can be produce by the @peterjc tools
[Galaxy wrappers for NCBI BLAST+ and related BLAST tools](https://github.com/peterjc/galaxy_blast)


Historic contributors
---------------------

 - Emma Prudent - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [UPMC](www.upmc.fr)/[CNRS](www.cnrs.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
 - Caroline Vernette - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [UPMC](www.upmc.fr)/[CNRS](www.cnrs.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
 - Gildas Le Corguillé @lecorguille - [ABiMS](http://abims.sb-roscoff.fr/) / [IFB](http://www.france-bioinformatique.fr/) - [UPMC](www.upmc.fr)/[CNRS](www.cnrs.fr) - [Station Biologique de Roscoff](http://www.sb-roscoff.fr/) - France
