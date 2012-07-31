<!-- This is the project specific website template -->
<!-- It can be changed as liked or replaced by other content -->
<?php

$domain=ereg_replace('[^\.]*\.(.*)$','\1',$_SERVER['HTTP_HOST']);
$group_name=ereg_replace('([^\.]*)\..*$','\1',$_SERVER['HTTP_HOST']);
$themeroot='https://r-forge.r-project.org/themes/rforge/';

echo '<?xml version="1.0" encoding="UTF-8"?>';
?>
<!DOCTYPE html
	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
	"http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en   ">
<head>
<meta http-equiv="Content-Type" content="text/html; charset=UTF-8" />
<title><?php echo $group_name; ?></title>
<link href="<?php echo $themeroot; ?>styles/estilo1.css" rel="stylesheet" type="text/css" />
<style type="text/css">
<!--
.style1 {font-size: small}
.R_code {
	font-family:"Courier New", Courier, monospace;
	font-size: x-small;
	font-style: italic;
}
.caption {
	font-size: x-small;
	font-style: italic;
	padding-bottom: 10px;
}
.R_env {
	font-family:"Courier New", Courier, monospace;
	color:#0000FF;
	font-size: x-small;
	margin-left: 10px;
}
.R_arg {
font-family:"Courier New", Courier, monospace;
color:#FF0000;
font-size: x-small;
}
-->
    </style>
</head>
<body>
<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
  <tr>
    <td><a href="http://r-forge.r-project.org/"><img src="<?php echo $themeroot; ?>imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td>
  </tr>
</table>
<!-- get project title  -->
<!-- own website starts here, the following may be changed as you like -->
<?php if ($handle=fopen('http://'.$domain.'/export/projtitl.php?group_name='.$group_name,'r')){
$contents = '';
while (!feof($handle)) {
	$contents .= fread($handle, 8192);
}
fclose($handle);
echo $contents; } ?>
<!-- end of project description -->
<div>
  <h1><strong>Tutorial: from soil profile data to 3D soil property and class maps</strong> (<a href="http://plotkml.r-forge.r-project.org/eberg.html">Eberg&ouml;tzen</a> case study)</h1>
</div>
<hr />
<p class="style1">Prepared by: <a href="http://www.wewur.wur.nl/popups/vcard.aspx?id=HENGL001" target="_blank">Tomislav Hengl</a>, <a href="http://www.wewur.wur.nl/popups/vcard.aspx?id=HEUVE015">Gerard B.M. Heuvelink</a> <br />
  Last update:
  <!-- #BeginDate format:Am1 -->July 31, 2012<!-- #EndDate -->
</p>
<p>The purpose of this tutorial is to demonstrate major processing steps used within the GSIF framework for generating soil property and soil class maps from point data, and with the help of multi-scale covariates. The GSIF (R package) <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. To learn more about the <strong>Global Soil Information Facilities</strong> (GSIF), visit the <a href="http://www.isric.org/projects/global-soil-information-facilities-gsif" target="_blank">main project page</a>. See also the complete list of <strong><a href="00Index.html">functions</a></strong> available via the GSIF package.</p>
<p>Download the tutorial as <a href="tutorial_eberg.R">R script</a>. </p>
<p><em><a name="top" id="top"></a>Table of content: </em></p>
<ul>
  <li><a href="#data_import">Loading the data and data screening</a>
    <ul>
      <li><a href="#histogram_plotting">Histogram plotting</a></li>
      <li><a href="#testing_CSR">Testing spatial randomness and point clustering</a></li>
      <li><a href="#MaxEnt_analysis">Testing the feature space coverage</a></li>
    </ul>
  </li>
  <li><a href="#model_fitting">Model fitting and spatial predictions</a>
    <ul>
      <li><a href="#preparing_point_data">Preparing the point data for geostatistical analysis</a></li>
      <li><a href="#generating_SPCs">Generating soil predictive components</a></li>
      <li><a href="#predicting_soil_properties">Predicting soil properties</a></li>
      <li><a href="#predicting_soil_classes">Predicting soil classes</a></li>
    </ul>
  </li>
  <li><a href="#multiscale_data">Predicting with multiscale data</a>
    <ul>
      <li><a href="#RK_multiscale">Regression-kriging using multiscale data</a></li>
      <li><a href="#RK_multisource">Merging multisource data</a></li>
    </ul>
  </li>
  <li><a href="#references">References</a></li>
</ul>
<hr />
<table width="100%" border="0" cellspacing="0" cellpadding="10">
  <tr>
    <th scope="col"><div align="left">
      <h2><a name="data_import" id="data_import"></a>Loading the data and  data screening</h2>
    </div></th>
    <th scope="col">&nbsp;</th>
    <th scope="col"><div align="right"><a href="#top">^to top</a> </div></th>
  </tr>
</table>
<p>For demonstration purposes we use the <a href="http://plotkml.r-forge.r-project.org/eberg.html">Eberg&ouml;tzen</a> case study, which   has been used by the SAGA GIS development team (<a href="http://saga-gis.org/en/about/references.html" target="_blank">Böhner et al., 2006</a>; <a href="http://dx.doi.org/10.1016/j.jag.2012.02.005" target="_blank">Hengl et al., 2012</a>). To start the exercise, first install and load all required packages (see also: <a href="http://gsif.r-forge.r-project.org">GSIF installation instructions</a>):</p>
<p class="R_code">&gt; library(plotKML)<br />
&gt; library(GSIF)</p>
<p class="R_env">Loading required package: RCurl<br />
Loading required package: bitops<br />
GSIF version 0.2-1 (2012-07-19)<br />
URL: http://gsif.r-forge.r-project.org/</p>
<p>load the input data:</p>
<p class="R_code">&gt; data(eberg)<br />
  &gt; data(eberg_grid)<br />
  &gt; data(eberg_grid25)</p>
<p>The <a href="http://plotkml.r-forge.r-project.org/eberg.html">Eberg&ouml;tzen</a> case study is a 10 by 10 km large case study in central Germany. The eberg object contains soil profile observations, and <span class="R_code">eberg_grid</span> and <span class="R_code">eberg_grid25</span> contains covariate layers at 100 and 25 m resolution. Before we proceed with running geostatistical analysis, we can run some initial data screening. We can first look at the data structure:</p>
<p class="R_code">&gt; str(eberg)</p>
<p class="R_env">'data.frame':   3670 obs. of  30 variables:<br />
  $ ID      : Factor w/ 3692 levels &quot;id0001&quot;,&quot;id0002&quot;,..: 3302 93 2827 1858 3539 3540 2828 94 1859 95 ...<br />
  $ soiltype: Factor w/ 13 levels &quot;A&quot;,&quot;B&quot;,&quot;D&quot;,&quot;G&quot;,..: NA NA NA NA NA NA NA NA NA NA ...<br />
  $ TAXGRSC : Factor w/ 13 levels &quot;Auenboden&quot;,&quot;Braunerde&quot;,..: NA NA NA NA NA NA NA NA NA NA ...<br />
  $ X       : int  3569323 3569328 3569328 3569335 3569336 3569340 3569343 3569344 3569348 3569349 ...<br />
  $ Y       : int  5716216 5715647 5716024 5715770 5716095 5716233 5716325 5715447 5714342 5716281 ...<br />
  $ UHDICM_A: num  0 0 0 0 0 0 0 0 0 0 ...<br />
  $ LHDICM_A: num  10 10 10 10 10 10 10 10 10 10 ...<br />
  $ SNDMHT_A: num  20 20 18.8 20 16.6 12.5 20 20 20 16.5 ...<br />
  $ SLTMHT_A: num  40 40 37.6 40 33.2 25 40 40 40 66.5 ...<br />
  $ CLYMHT_A: num  40 40 37.6 40 33.2 25 40 40 40 21 ...<br />
  $ UHDICM_B: num  10 10 10 10 10 10 10 10 10 10 ...<br />
  $ LHDICM_B: num  30 30 30 30 30 30 30 30 30 30 ...<br />
  $ SNDMHT_B: num  NA 20 NA NA NA NA NA 20 18.8 16.5 ...<br />
  $ SLTMHT_B: num  NA 40 NA NA NA NA NA 40 37.6 66.5 ...<br />
  $ CLYMHT_B: num  NA 40 NA NA NA NA NA 40 37.6 21 ...<br />
  $ UHDICM_C: num  30 30 30 30 30 30 30 30 30 30 ...<br />
  $ LHDICM_C: num  50 50 50 50 50 50 50 50 50 50 ...<br />
  $ SNDMHT_C: num  NA 18.8 NA NA NA NA NA 18.8 NA NA ...<br />
  $ SLTMHT_C: num  NA 37.6 NA NA NA NA NA 37.6 NA NA ...<br />
  $ CLYMHT_C: num  NA 37.6 NA NA NA NA NA 37.6 NA NA ...<br />
  $ UHDICM_D: num  50 50 50 50 50 50 50 50 50 50 ...<br />
  $ LHDICM_D: num  70 70 70 70 70 70 70 70 70 70 ...<br />
  $ SNDMHT_D: num  NA 18.8 NA NA NA NA NA 18.8 NA NA ...<br />
  $ SLTMHT_D: num  NA 37.6 NA NA NA NA NA 37.6 NA NA ...<br />
  $ CLYMHT_D: num  NA 37.6 NA NA NA NA NA 37.6 NA NA ...<br />
  $ UHDICM_E: num  70 70 70 70 70 70 70 70 70 70 ...<br />
  $ LHDICM_E: num  90 90 90 90 90 90 90 90 90 90 ...<br />
  $ SNDMHT_E: num  NA 18.8 NA NA NA NA NA 18.8 NA NA ...<br />
  $ SLTMHT_E: num  NA 37.6 NA NA NA NA NA 37.6 NA NA ...<br />
$ CLYMHT_E: num  NA 37.6 NA NA NA NA NA 37.6 NA NA ... </p>
<p>This shows that soil texture fractions <span class="R_code">SNDMHT</span>, <span class="R_code">SLTMHT</span> and <span class="R_code">CLYMHT</span> have been sampled at fixed depths (10-30, 30-50, 50-70, 70-90), and <span class="R_code">soiltype</span> contains <a href="http://plotkml.r-forge.r-project.org/eberg.html">13 classes of soil types</a>. </p>
<h3><a name="histogram_plotting" id="histogram_plotting"></a>Histogram plotting </h3>
<p>Our focus in this exercises is spatial prediction of sand content over the whole volume of soil. Let us see some summary info for sand content:</p>
<p class="R_code">&gt; summary(eberg$SNDMHT_A)</p>
<p class="R_env">Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's <br />
  4.062  18.800  20.000  30.260  42.300  92.500   3.000 </p>
<p class="R_code">&gt; library(StatDA)<br />
  &gt; par(mar=c(2.5,2.5,0.5,0.5), oma=c(0,0,0,0))<br />
  &gt; edaplot(eberg$SNDMHT_A[!is.na(eberg$SNDMHT_A)], H.freq=TRUE, box=FALSE, S.pch=3, S.cex=0.5, D.lwd=1.5, P.ylab=&quot;&quot;, P.log=FALSE, P.logfine=c(5,10), P.main=&quot;&quot;, P.xlab=&quot;&quot;, B.pch=3, B.cex=0.5)</p>
<table width="350" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Histogram for sand content.
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_hist_SNDMHT.png" alt="Fig_eberg_hist_SNDMHT.png" width="350" /></th>
  </tr>
</table>
<p>From the plot above, we can also observe that many numbers are in fact overlapping. It seems that there is only a clusters of 
  values possible for <span class="R_code">SNDMHT</span>. Going 
  back to the origin of this data set, we can notice that the sand, silt and clay values have 
  been determined by using the so-called texture by hand method, i.e. via texture classes. The literature (<a href="https://www.soils.org/publications/sssaj/articles/65/4/1038" target="_blank">Skaggs
  et al., 2001</a>) shows that the this technique can be used to determine the content of soil
earth fractions only to an accuracy of ±5–10%. This means that we should not plan to map any of the texture fractions to a precision better than ±5% (detection limit), because it would exceed the measurement error.</p>
<p>These is a relatively large set (3670 points), so it might be a good idea to subset it (e.g. take only 30% of samples) to speed up the processing:</p>
<p class="R_code">&gt; eberg.xy &lt;- eberg[runif(nrow(eberg)) &lt; .3,]<br />
  &gt; coordinates(eberg.xy) &lt;- ~X+Y<br />
  &gt; proj4string(eberg.xy) &lt;- CRS(&quot;+init=epsg:31467&quot;)<br />
  &gt; gridded(eberg_grid) &lt;- ~x+y<br />
  &gt; proj4string(eberg_grid) &lt;- CRS(&quot;+init=epsg:31467&quot;)<br />
  &gt; gridded(eberg_grid25) &lt;- ~x+y<br />
  &gt; proj4string(eberg_grid25) &lt;- CRS(&quot;+init=epsg:31467&quot;)</p>
<h3><a name="testing_CSR" id="testing_CSR"></a>Testing spatial randomness and point clustering</h3>
<p>Prior to geostatistical analysis, it is usefull to assess how representative and consistent 
  the soil profile data is. To achieve this, we can run some basic analysis of the point geometry and 
  then overlap the points with predictors to see how well are the environmental features 
  represented. First, we can test if the points represent the geographical space in an unbiased way (see chapter 7, in <a href="http://asdar-book.org">Bivand et al., 2008</a>):</p>
<p class="R_code">&gt; library(spatstat) <br />
&gt; eberg.ppp &lt;- ppp(x=coordinates(eberg.xy)[,1], y=coordinates(eberg.xy)[,2], marks=eberg.xy$zinc, window=mg_owin)<br />
&gt; env.eberg.xy &lt;- envelope(eberg.ppp, fun=Gest)</p>
<p class="R_env">Min. 1st Qu.  Median    Mean 3rd Qu.    Max. <br />
  17.80   78.22  120.70  140.00  186.80  741.10 &gt; env.eberg.xy &lt;- envelope(eberg.ppp, fun=Gest)<br />
  Generating 99 simulations of CSR  ...<br />
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,<br />
  16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,<br />
  31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,<br />
  46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,<br />
  61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75,<br />
  76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,<br />
  91, 92, 93, 94, 95, 96, 97, 98, 99.<br />
Done.</p>
<p class="R_code">&gt; par(mar=c(4.5,4.5,0.5,0.5), oma=c(0,0,0,0))<br />
&gt; plot(env.eberg.xy, lwd=list(3,1,1,1), main=&quot;&quot;)</p>
<p class="R_env">lty col  key      label                                           meaning<br />
  obs    1   1  obs  G[obs](r)           observed value of G(r) for data pattern<br />
  theo   2   2 theo G[theo](r)                 theoretical value of G(r) for CSR<br />
  hi     1   8   hi   G[hi](r) upper pointwise envelope of G(r) from simulations<br />
  lo     1   8   lo   G[lo](r) lower pointwise envelope of G(r) from simulations</p>
<table width="350" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: The results of the Complete Spatial Randomness test (spatstat package) for the Eberg&ouml;tzen case study.
  In this case the observed G-curve does not fit into the 95% range for the CSR, for this study area. 
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_CRS_test.png" alt="Fig_eberg_CRS_test.png" width="350" /></th>
  </tr>
</table>
<p>which shows that the point samples do not exactly satisfy the Complete Spatial Randomness test (<a href="http://cran.r-project.org/web/packages/spatstat/spatstat.pdf" target="_blank">spatstat</a> package) &#8212; samples at larger distances have a lower spatial density than a completely random design, which usually means that some parts of the case study are under-sampled (see also figure below). On the other hand, points at shorter distances (&lt;75 m) would pass the CSR test, which indicates that, at least at shorter distances, there is no clustering in geographical space. </p>
<h3><a name="MaxEnt_analysis" id="MaxEnt_analysis"></a>Testing the feature space coverage</h3>
<p>To see how representative is the point data considering the coverage of feature space, we use the <a href="MaxEnt-method.html">MaxEnt</a> function that extends methods available via the <a href="http://cran.r-project.org/web/packages/dismo/">dismo</a> package:</p>
<p class="R_code">&gt; me.eberg &lt;- MaxEnt(occurrences=eberg.ppp, covariates=eberg_grid)<br />
  &gt; par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,0.5), oma=c(0,0,0,0))<br />
  &gt; image(as(me.eberg@predicted, &quot;SpatialPixelsDataFrame&quot;), col=rev(heat.colors(25)), xlab=&quot;&quot;, ylab=&quot;&quot;)<br />
  &gt; points(me.eberg@occurrences, pch=&quot;+&quot;, cex=.7)<br />
  &gt; image(me.eberg@sp.domain, col=&quot;grey&quot;, xlab=&quot;&quot;, ylab=&quot;&quot;)</p>
<table width="500" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Feature space analysis: sampling likelihood in feature space based on MaxEnt (left; dark red indicates higher values), and areas fully represented by the current samples following the cross-validation (right). See code examples for more details.
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_MaxEnt_test.png" alt="Fig_eberg_MaxEnt_test.png" width="500" /></th>
  </tr>
</table>
<p>Which shows that, some parts of the study area (higher elevations, specific land cover types) have been systematically omitted from sampling. The map on the right shows which pixels (grey) are actually valid to run predictions as the probability of occurrence of sampling points at these locations (at least based on the feature space defined by the covariates) is statistically significant. Ideally, the whole map should be gray, which means none of the areas in feature space have been under-represented. This is nothing that should worry us too much, but something we need to be aware when doing the interpretation of produced maps. To learn  more about <a href="http://www.cs.princeton.edu/~schapire/maxent/" target="_blank">MaxEnt</a>, refer to the dismo package vignette (<a href="http://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf" target="_blank">Hijmans and Elith, 2012</a>). </p>
<hr />
<table width="100%" border="0" cellpadding="10" cellspacing="0">
  <tr>
    <th scope="col"><div align="left">
      <h2><a name="model_fitting" id="model_fitting"></a>Model fitting and spatial predictions</h2>
    </div></th>
    <th scope="col">&nbsp;</th>
    <th scope="col"><div align="right"><a href="#top">^to top</a> </div></th>
  </tr>
</table>
<p>Any geostatistical mapping boils down to two major processes: (<em>a</em>) model fitting (or model estimation), and (<em>b</em>) generation of spatial predictions. Also in GSIF package geostatistical mapping is implemented via two main functions: <a href="fit.gstatModel-method.html">fit.gstatModel</a> and <a href="gstatModel-class.html">predict</a>. In the following sections we demonstrate how to fit models and generate predictions using soil profile data and a diversity of covariate layers typical for soil mapping applications. Functions <a href="fit.gstatModel-method.html">fit.gstatModel</a> and <a href="gstatModel-class.html">predict</a> largely extend the functionality of the <a href="http://stat.ethz.ch/R-manual/R-patched/library/stats/html/00Index.html">stats</a> and <a href="http://www.gstat.org">gstat</a> packages. Gstat package (see chapter 8, in <a href="http://asdar-book.org">Bivand et al., 2008</a>) is by many still considered one of the most comprehensive packages for geostatistical modeling.</p>
<h3><a name="preparing_point_data" id="preparing_point_data"></a>Preparing point data for geostatistical analysis</h3>
<p>Before we build a regression-kriging model,  we need to prepare the soil observations to a format suitable for fitting of 3D models. We start by converting the <span class="R_code">eberg</span> data frame to class <a href="http://casoilresource.lawr.ucdavis.edu/drupal/node/966">SoilProfileCollection-class</a> (<a href="http://cran.r-project.org/package=aqp">aqp</a> package) and the <a href="geosamples-class.html">geosamples-class</a> (GSIF package):</p>
<p class="R_code">&gt; s.lst &lt;- c(&quot;ID&quot;, &quot;soiltype&quot;, &quot;TAXGRSC&quot;, &quot;X&quot;, &quot;Y&quot;)<br />
&gt; h.lst &lt;- c(&quot;UHDICM&quot;,&quot;LHDICM&quot;,&quot;SNDMHT&quot;,&quot;SLTMHT&quot;,&quot;CLYMHT&quot;)<br />
&gt; sites &lt;- eberg[,s.lst]<br />
&gt; horizons &lt;- getHorizons(eberg, idcol=&quot;ID&quot;, sel=h.lst)<br />
&gt; eberg.spc &lt;- join(horizons, sites, type='inner')</p>
<p class="R_env"> Joining by: ID</p>
<p class="R_code">&gt; depths(eberg.spc) &lt;- ID ~ UHDICM + LHDICM<br />
&gt; site(eberg.spc) &lt;- as.formula(paste(&quot;~&quot;, paste(s.lst[-1], collapse=&quot;+&quot;), sep=&quot;&quot;))</p>
<p class="R_env">Warning message:<br />
converting IDs from factor to character</p>
<p class="R_code">&gt; coordinates(eberg.spc) &lt;- ~X+Y<br />
&gt; proj4string(eberg.spc) &lt;- CRS(&quot;+init=epsg:31467&quot;)<br />
&gt; # convert to logits: <br />
&gt; eberg.spc@horizons$SNDMHT.t &lt;- log((eberg.spc@horizons$SNDMHT/100)/(1-eberg.spc@horizons$SNDMHT/100))</p>
<p>where <span class="R_code">SNDMHT.t</span> is the logit-transformed value of the target variable. Logit transformation is required to prevent from making predictions outside the 0-1 range. </p>
<p>GSIF package by default works with <a href="geosamples-class.html">geosamples-class</a> as the main class for field observations (points). Conversion of the SPC class data to geosamples is straight forward:</p>
<p class="R_code">&gt; eberg.geo &lt;- as.geosamples(eberg.spc)</p>
<p class="R_env">  Reprojecting to +proj=longlat +datum=WGS84 ...</p>
<p class="R_code">&gt; str(eberg.geo)</p>
<p class="R_env">Formal class 'geosamples' [package &quot;GSIF&quot;] with 3 slots<br />
..@ registry: chr NA<br />
..@ methods :'data.frame':    6 obs. of  4 variables:<br />
.. ..$ methodid      : Factor w/ 6 levels &quot;CLYMHT&quot;,&quot;SLTMHT&quot;,..: 1 2 3 4 5 6<br />
.. ..$ description   : logi [1:6] NA NA NA NA NA NA<br />
.. ..$ units         : logi [1:6] NA NA NA NA NA NA<br />
.. ..$ detectionLimit: logi [1:6] NA NA NA NA NA NA<br />
..@ data    :'data.frame':    25762 obs. of  14 variables:<br />
.. ..$ observationid   : chr [1:25762] NA NA NA NA ...<br />
.. ..$ sampleid        : chr [1:25762] &quot;id0003&quot; &quot;id0007&quot; &quot;id0014&quot; &quot;id0015&quot; ...<br />
.. ..$ longitude       : num [1:25762] 10.2 10.1 10.1 10.1 10.1 ...<br />
.. ..$ latitude        : num [1:25762] 51.6 51.5 51.5 51.6 51.6 ...<br />
.. ..$ locationError   : num [1:25762] NA NA NA NA NA NA NA NA NA NA ...<br />
.. ..$ TimeSpan.begin  : POSIXct[1:25762], format: NA NA NA ...<br />
.. ..$ TimeSpan.end    : POSIXct[1:25762], format: NA NA NA ...<br />
.. ..$ altitude        : num [1:25762] 0 0 0 0 0 0 0 0 0 0 ...<br />
.. ..$ altitudeMode    : chr [1:25762] &quot;relativeToGround&quot; &quot;relativeToGround&quot; &quot;relativeToGround&quot; &quot;relativeToGround&quot; ...<br />
.. ..$ sampleArea      : num [1:25762] 1 1 1 1 1 1 1 1 1 1 ...<br />
.. ..$ sampleThickness : num [1:25762] 2 2 2 2 2 2 2 2 2 2 ...<br />
.. ..$ observedValue   : chr [1:25762] &quot;L&quot; &quot;L&quot; &quot;B&quot; &quot;K&quot; ...<br />
.. ..$ methodid        : Factor w/ 6 levels &quot;CLYMHT&quot;,&quot;SLTMHT&quot;,..: 5 5 5 5 5 5 5 5 5 5 ...<br />
.. ..$ measurementError: num [1:25762] NA NA NA NA NA NA NA NA NA NA ...</p>
<p class="R_code">&gt; levels(eberg.geo@data$methodid)</p>
<p class="R_env">  [1] &quot;CLYMHT&quot;   &quot;SLTMHT&quot;   &quot;SNDMHT&quot;   &quot;SNDMHT.t&quot; &quot;soiltype&quot;<br />
  [6] &quot;TAXGRSC&quot;</p>
<p><a href="geosamples-class.html">Geosamples-class</a> can be considered the most plain (standard) format for any space-time observations and, as such, highly suitable for storing free-form geosamples. Main advantage of using this class in R is that it can be easily manipulated and converted to spatial and aqp classes. Likewise, the column names in the <span class="R_code">@data</span> slot correspond to the tag names used in the KML schema, which makes it easier to export such data to some GIS or Google Earth. Note also that the idea of using <a href="geosamples-class.html">geosamples</a> is to allow users to extend possibilities of geostatistical analysis with data parameters such as horizontal and vertical support (<span class="R_code">sampleArea</span> and <span class="R_code">sampleThikness</span>) and uncertainty measures (<span class="R_code">locationError</span> and <span class="R_code">measurementError</span>).</p>
<h3><a name="generating_SPCs" id="generating_SPCs"></a>Generating soil predictive components</h3>
<p>Prior to geostatistical modeling, it is also probably a good idea to convert all covariates to independent components. This way, it will be easier to subset to the optimal number of predictors during the regression analysis. PCA helps reducing the prediction bias, which might happen if the covariates are cross-correlated. A wrapper function <a href="spc-method.html">spc</a> will convert all factor variables to indicators and run PCA on a stack of grids: </p>
<p class="R_code">&gt; formulaString &lt;- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6<br />
&gt; eberg_spc &lt;- spc(eberg_grid, formulaString)</p>
<p class="R_env">Converting PRMGEO6 to indicators...<br />
  Converting covariates to principal components...</p>
<p>The output is a stack of independent components, all numeric and all scaled around 0 value. To see which inputs define any of the components, we can look at the rotation table or plot the images: </p>
<p class="R_code">&gt; eberg_spc@pca$rotation<br />
&gt; pal = rev(rainbow(65)[1:48])<br />
  &gt; rd = range(eberg_spc@predicted@data[,1], na.rm=TRUE)<br />
&gt; spplot(eberg_spc@predicted[1:4], at=seq(rd[1], rd[2], length.out=48), col.regions=pal)</p>
<table width="500" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: First four components derived using the eberg_grid data.
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_SPCs1_4.png" alt="Fig_eberg_SPCs1_4.png" width="500" /></th>
  </tr>
</table>
<h3><a name="predicting_soil_properties" id="predicting_soil_properties"></a>Predicting soil properties</h3>
<p>In the GSIF package, regression-kriging model can be fitted at once by using the <a href="fit.gstatModel-method.html">fit.gstatModel</a> function. First, we need to specify the model explaining the distribution of soil property as a function of soil covariates: </p>
<p class="R_code">&gt; glm.formulaString = as.formula(paste(&quot;observedValue ~ &quot;, paste(names(eberg_spc@predicted), collapse=&quot;+&quot;), &quot;+ ns(altitude, df=4)&quot;))<br />
&gt; glm.formulaString</p>
<p class="R_env">  observedValue ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + <br />
PC9 + PC10 + PC11 + ns(altitude, df = 4)</p>
<p>In other words the observed values will be modeled as a function of PCs and altitude (natural splines via the <a href="http://stat.ethz.ch/R-manual/R-devel/library/splines/html/ns.html">ns</a> function). In the GSIF package, the 3D GLM-kriging model can be fitted at once by running: </p>
<p class="R_code">&gt; SNDMHT.m &lt;- fit.gstatModel(observations=eberg.geo, glm.formulaString, covariates=eberg_spc@predicted, methodid=&quot;SNDMHT.t&quot;)<br />
 &gt; summary(SNDMHT.m@regModel)</p>
<p class="R_env">... </p>
<p class="R_env">Coefficients:<br />
  Estimate Std. Error t value<br />
  (Intercept)           -1.231e+00  2.119e-02 -58.086<br />
  PC1                   -4.030e-01  2.247e-02 -17.937<br />
  PC2                   -1.741e-01  9.001e-03 -19.343<br />
  PC4                    4.909e-01  6.209e-03  79.061<br />
  PC5                    1.262e-01  6.260e-03  20.154<br />
  PC6                    3.566e-02  1.055e-02   3.381<br />
  PC7                    6.567e-02  9.102e-03   7.214<br />
  PC8                    4.918e-01  1.988e-02  24.745<br />
  PC9                   -3.532e-01  2.202e-02 -16.042<br />
  PC10                   5.050e-01  3.535e-02  14.284<br />
  PC11                  -1.809e+14  1.202e+13 -15.044<br />
  ns(altitude, df = 4)1  5.857e-02  3.448e-02   1.698<br />
  ns(altitude, df = 4)2  1.475e-01  2.894e-02   5.097<br />
  ns(altitude, df = 4)3  2.406e-01  4.622e-02   5.205<br />
  ns(altitude, df = 4)4  2.270e-01  2.329e-02   9.743<br />
  Pr(&gt;|t|) <br />
  (Intercept)            &lt; 2e-16 ***<br />
  PC1                    &lt; 2e-16 ***<br />
  PC2                    &lt; 2e-16 ***<br />
  PC4                    &lt; 2e-16 ***<br />
  PC5                    &lt; 2e-16 ***<br />
  PC6                   0.000724 ***<br />
  PC7                   5.75e-13 ***<br />
  PC8                    &lt; 2e-16 ***<br />
  PC9                    &lt; 2e-16 ***<br />
  PC10                   &lt; 2e-16 ***<br />
  PC11                   &lt; 2e-16 ***<br />
  ns(altitude, df = 4)1 0.089448 . <br />
  ns(altitude, df = 4)2 3.50e-07 ***<br />
  ns(altitude, df = 4)3 1.97e-07 ***<br />
  ns(altitude, df = 4)4  &lt; 2e-16 ***<br />
  ---<br />
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 </p>
<p class="R_env">(Dispersion parameter for gaussian family taken to be 0.6365825)</p>
<p class="R_env"> Null deviance: 13537.1  on 11922  degrees of freedom<br />
  Residual deviance:  7580.4  on 11908  degrees of freedom<br />
  (1967 observations deleted due to missingness)<br />
  AIC: 28468</p>
<p class="R_env">Number of Fisher Scoring iterations: 2</p>
<p class="R_code">&gt; SNDMHT.m@vgmModel</p>
<p class="R_env">  model     psill    range kappa ang1 ang2 ang3 anis1<br />
  1   Nug 0.2771894   0.0000   0.0    0    0    0     1<br />
  2   Exp 0.4396340 145.4435   0.5    0    0    0     1<br />
  anis2<br />
  1 1.00000<br />
  2 0.00015</p>
<p>These result show that the model is significant, and this is valid for both GLM and the variogram. Note however that this is not completely a 3D regression-krging, as we do not actually have values of PCs at different depths (in fact, most of PCs relate only to the surface), so that many values of covariates are basically copied to all depths. This does not represent any problem for the GLM modeling, however, you should be aware that, because values of covariates are fixed with the different depths, the resulting 3D patterns in the target variable will be mainly controlled by the surface patterns.</p>
<p>Once we have fitted a <a href="gstatModel-class.html">gstatModel</a>, we can generate predictions and estimate the associated uncertainty at any depth. In the last step, we need to prepare the 3D prediction locations i.e. grid cells that need to be mapped. In GSIF package, this can be done by using the <a href="make.3Dgrid-method.html">sp3D</a> function:</p>
<p class="R_code">&gt; new3D &lt;- sp3D(eberg_spc@predicted)</p>
<p>This will prepare the existing SPCs for 3D predictions (6 standard depths, at block support responding to the size of standard horizons) by adding the <span class="R_code">altitude</span> column:</p>
<p class="R_code">&gt; str(new3D[[1]]@grid)</p>
<p class="R_env">  Formal class 'GridTopology' [package &quot;sp&quot;] with 3 slots<br />
  ..@ cellcentre.offset: Named num [1:3] 3.57e+06 5.71e+06 -2.50e-02<br />
  .. ..- attr(*, &quot;names&quot;)= chr [1:3] &quot;longitude&quot; &quot;latitude&quot; &quot;altitude&quot;<br />
  ..@ cellsize         : Named num [1:3] 1e+02 1e+02 5e-02<br />
  .. ..- attr(*, &quot;names&quot;)= chr [1:3] &quot;longitude&quot; &quot;latitude&quot; &quot;altitude&quot;<br />
  ..@ cells.dim        : Named int [1:3] 100 100 1<br />
  .. ..- attr(*, &quot;names&quot;)= chr [1:3] &quot;longitude&quot; &quot;latitude&quot; &quot;altitude&quot; </p>
<p>So that we can now run predictions for each standard depth in a loop: </p>
<p class="R_code">&gt; sd.l &lt;- lapply(new3D, FUN=function(x){predict(SNDMHT.m, predictionLocations=x, nfold=0)})</p>
<p class="R_env">Generating predictions using the trend model (KED method)...<br />
  [using universal kriging]</p>
<p class="R_env">...</p>
<p>This operation can take time depending on the size of the grids and number of 3D points used to generate predictions. </p>
<p>Finally, we can prepare the produced predictions and export them as <a href="GlobalSoilMap-class.html">GlobalSoilMap-class</a> object. First, we back-transform the predictions to the 0-100% scale:</p>
<p class="R_code">&gt; for(j in 1:length(sd.l)){ sd.l[[j]]@predicted$observedValue &lt;- exp(sd.l[[j]]@predicted$observedValue)/(1+exp(sd.l[[j]]@predicted$observedValue))*100 }</p>
<p>and then  reproject the produced predictions to  geographical coordinates (<a href="http://spatialreference.org/ref/epsg/4326/" target="_blank">WGS84</a>) using the <a href="make.3Dgrid-method.html">make.3Dgrid</a> function: </p>
<p class="R_code">&gt; p = get(&quot;cellsize&quot;, envir = GSIF.opts)[2]<br />
  &gt; s = get(&quot;stdepths&quot;, envir = GSIF.opts)<br />
  &gt; s</p>
<p class="R_env">  [1] -0.025 -0.075 -0.225 -0.450 -0.800 -1.500</p>
<p class="R_code">&gt; sd.ll &lt;- sapply(1:length(sd.l), FUN=function(x){make.3Dgrid(sd.l[[x]]@predicted[3:4], pixelsize=p, stdepths=s[x])})</p>
<p class="R_env">  Resampling 2 layers to +proj=longlat +datum=WGS84 with grid cell size of 0.000833333333333333 ...</p>
<p>The object <span class="R_code">sd.ll</span> can now be saved as a <a href="GlobalSoilMap-class.html">GlobalSoilMap-class</a> object:</p>
<p class="R_code">&gt; SNDMHT.gsm &lt;- GlobalSoilMap(varname=&quot;SNDMHT&quot;, sd.ll, period=c(&quot;1999-02-01&quot;, &quot;2001-07-01&quot;))<br />
&gt; str(SNDMHT.gsm, max.level=2)</p>
<p class="R_env">Formal class 'GlobalSoilMap' [package &quot;.GlobalEnv&quot;] with 9 slots<br />
  ..@ varname       : chr &quot;SNDMHT&quot;<br />
  ..@ TimeSpan.begin: POSIXct[1:1], format: &quot;1999-02-01&quot;<br />
  ..@ TimeSpan.end  : POSIXct[1:1], format: &quot;2001-07-01&quot;<br />
  ..@ sd1           :Formal class 'SpatialPixelsDataFrame' [package &quot;sp&quot;] with 7 slots<br />
  ..@ sd2           :Formal class 'SpatialPixelsDataFrame' [package &quot;sp&quot;] with 7 slots<br />
  ..@ sd3           :Formal class 'SpatialPixelsDataFrame' [package &quot;sp&quot;] with 7 slots<br />
  ..@ sd4           :Formal class 'SpatialPixelsDataFrame' [package &quot;sp&quot;] with 7 slots<br />
  ..@ sd5           :Formal class 'SpatialPixelsDataFrame' [package &quot;sp&quot;] with 7 slots<br />
..@ sd6           :Formal class 'SpatialPixelsDataFrame' [package &quot;sp&quot;] with 7 slots</p>
<p class="R_code"><span class="R_code">&gt; save(SNDMHT.gsm, file=&quot;SNDMHT.rda&quot;, compress=&quot;xz&quot;)</span>  </p>
<p>and visualized  in Google Earth by using the <a href="http://plotkml.r-forge.r-project.org/">plotKML</a> package functionality: </p>
<p class="R_code">&gt;  z0 = mean(eberg_grid$DEMSRT6, na.rm=TRUE)<br />
  &gt; for(j in 1:length(sd.ll)){<br />
  + kml(slot(SNDMHT.gsm, paste(&quot;sd&quot;, j, sep=&quot;&quot;)), folder.name = paste(&quot;eberg_sd&quot;, j, sep=&quot;&quot;), <br />
  + 
  file = paste(&quot;SNDMHT_sd&quot;, j, &quot;.kml&quot;, sep=&quot;&quot;), colour = observedValue, zlim=c(10,85), <br />
  + 
  raster_name = paste(&quot;SNDMHT_sd&quot;, j, &quot;.png&quot;, sep=&quot;&quot;), altitude = z0+5000+(s[j]*2500))<br />
  + }</p>
<p class="R_env">  KML file header opened for parsing...<br />
  Parsing to KML...<br />
  Closing  SNDMHT_sd1.kml<br />
  KML file header opened for parsing...<br />
  Parsing to KML...<br />
  Closing  SNDMHT_sd2.kml<br />
  KML file header opened for parsing...<br />
  Parsing to KML...<br />
  Closing  SNDMHT_sd3.kml<br />
  KML file header opened for parsing...<br />
  Parsing to KML...<br />
  Closing  SNDMHT_sd4.kml<br />
  KML file header opened for parsing...<br />
  Parsing to KML...<br />
  Closing  SNDMHT_sd5.kml<br />
  KML file header opened for parsing...<br />
  Parsing to KML...<br />
  Closing  SNDMHT_sd6.kml</p>
<table width="500" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Predicted sand content using 3D GLM-kriging as visualized in Google Earth (<a href="eberg_SNDMHT_6depths.zip">zip</a>). 
  When visualizing this data, make sure that the relief exaggeration in Google Earth has been set at 1. 
  </caption>
  <tr>
    <th scope="col"><a href="eberg_SNDMHT_6depths.zip"><img src="Fig_GSIF_predicted_layers.jpg" alt="Fig_GSIF_predicted_layers.jpg" width="500" border="0" /></a></th>
  </tr>
</table>
<h3><a name="predicting_soil_classes" id="predicting_soil_classes"></a>Predicting soil classes</h3>
<p>GSIF package also provides functionality for pedometric mapping of soil classes. Soil types can be mapped using a wrapper function   <a href="spfkm.html">spfkm</a>. This will   run supervised fuzzy <em>k</em>-means using a list of covariates layers provided as <span class="R_code">&quot;SpatialPixelsDataFrame&quot;</span> object, and optional class centers and class variances. As in the case of continuous/numeric variables, the process consists of model fitting and predictions. In GSIF package, the two are wrapped into a single function: </p>
<p class="R_code">&gt; eberg_sm &lt;- spfkm(formulaString, eberg.xy, eberg_spc@predicted)</p>
  <span class="R_env">Loading required package: nnet<br />
Fitting a multinomial logistic regression model...<br />
# weights:  144 (121 variable)<br />
initial  value 1856.225267 <br />
iter  10 value 1196.619912<br />
iter  20 value 1166.971099<br />
iter  30 value 1156.198652<br />
iter  40 value 1150.212361<br />
iter  50 value 1145.896491<br />
iter  60 value 1142.037781<br />
iter  70 value 1141.479941<br />
iter  80 value 1140.839636<br />
iter  90 value 1140.699639<br />
final  value 1140.698143 <br />
converged<br />
Loading required package: mda<br />
Loading required package: class</span>
<p class="R_env">Attaching package: ‘class’</p>
<p class="R_env">The following object(s) are masked from ‘package:reshape’:</p>
<p class="R_env"> condense</p>
<p class="R_env">Estimated prediction error: 0.5932</p>
<p>The output is an object of class <a href="SpatialMemberships-class.html">SpatialMemberships</a>, which contains predictions, the model parameters, class centres and variances, and membership maps derived per class. In the case above,  the class centers and variances were not specified, hence <a href="spfkm.html">spfkm</a> tries to estimate them by fitting a multinomial logistic regression model (<a href="http://stat.ethz.ch/R-manual/R-devel/library/nnet/html/multinom.html" target="_blank">multinom</a>) available via the package <a href="http://cran.r-project.org/package=nnet" target="_blank">nnet</a> (<a href="http://www.stats.ox.ac.uk/pub/MASS4/">Venables and Ripley, 2002</a>). Class centres and variances can also be passed via the <span class="R_code">class.c</span> and <span class="R_code">class.sd</span> arguments. in which case <a href="spfkm.html">spfkm</a> will directly derived supervised fuzzy <em>k</em>-means. To see the estimated class centers use e.g.:</p>
<p class="R_code">&gt; eberg_sm@class.c[1,]</p>
<p class="R_env"> PC1        PC2        PC3        PC4        PC5        PC6        PC7        PC8 <br />
3.1607835  1.7982977  0.2260035 -0.5133303 -0.9440280 -0.2689071  0.4962124  1.9947157 <br />
PC9       PC10 <br />
0.8995294 -0.3908947</p>
<p class="R_code">&gt; row.names(eberg_sm@class.c) </p>
<p class="R_env">[1] &quot;A&quot;  &quot;B&quot;  &quot;D&quot;  &quot;G&quot;  &quot;Hw&quot; &quot;L&quot;  &quot;N&quot;  &quot;Q&quot;  &quot;R&quot;  &quot;S&quot;  &quot;Z&quot;</p>
<table width="500" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Predicted soil types for the Eberg&ouml;tzen case study.
  See <a href="spfkm.html">spfkm</a> for more details. This can be also considered a downscaling or disaggregation method. 
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_Soiltypes_spfkm.png" alt="Fig_eberg_Soiltypes_spfkm.png" width="500" /></th>
  </tr>
</table>
<p>Mapped soil class memberships can also be used to map soil properties. The regression model changes to e.g.:</p>
<p class="R_code">&gt; glm.formulaString2</p>
<p class="R_env">  SNDMHT_A ~ A + B + D + G + Hw + L + N + Q + R + S + Z - 1</p>
<p class="R_code">&gt; SNDMHT.m2 &lt;- fit.gstatModel(observations=eberg.xy, glm.formulaString2, covariates=eberg_sm@mu)<br />
 &gt; summary(SNDMHT.m2@regModel)</p>
<p class="R_env">Call:<br />
  glm(formula = SNDMHT_A ~ A + B + D + Hw + L + N + Q + R + S + <br />
Z - 1, family = family, data = x)</p>
<p class="R_env">Deviance Residuals: <br />
  Min       1Q   Median       3Q      Max <br />
  -39.449   -9.981   -1.528    8.040   65.483 </p>
<p class="R_env">Coefficients:<br />
  Estimate Std. Error t value Pr(&gt;|t|) <br />
  A    17.194     10.972   1.567    0.117 <br />
  B    50.183      1.399  35.881  &lt; 2e-16 ***<br />
  D    12.473      2.660   4.689 3.22e-06 ***<br />
  Hw   16.842      2.871   5.866 6.47e-09 ***<br />
  L    25.997      1.244  20.898  &lt; 2e-16 ***<br />
  N    56.134      3.142  17.867  &lt; 2e-16 ***<br />
  Q    37.881      2.580  14.683  &lt; 2e-16 ***<br />
  R    12.037      8.429   1.428    0.154 <br />
  S    28.960      1.597  18.136  &lt; 2e-16 ***<br />
  Z    17.904      1.864   9.604  &lt; 2e-16 ***<br />
  ---<br />
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 </p>
<p class="R_env">(Dispersion parameter for gaussian family taken to be 292.8399)</p>
<p class="R_env"> Null deviance: 1191564  on 828  degrees of freedom<br />
  Residual deviance:  239543  on 818  degrees of freedom<br />
  (1 observation deleted due to missingness)<br />
  AIC: 7064.4</p>
<p class="R_env">Number of Fisher Scoring iterations: 2</p>
<p>Note that intercept needs to be taken out, so that the best predictor of the sand content for some soil type is basically the mean value of the sand for that soil type (for example class B is expected to have an average sand content of about 50.1%). If you compare the model based on soil classes and model fitted in the previous section (<span class="R_code">SNDMHT.m</span>), you can see that fitting data using a 3D model results in a slightly better fit. Nevertheless, soil classes are in this case study significant estimators of sand content.</p>
<hr />
<table width="100%" border="0" cellspacing="0" cellpadding="10">
  <tr>
    <th scope="col"><div align="left">
      <h2><a name="multiscale_data" id="multiscale_data"></a>Predicting with multi-scale data</h2>
    </div></th>
    <th scope="col">&nbsp;</th>
    <th scope="col"><div align="right"><a href="#top">^to top</a> </div></th>
  </tr>
</table>
<p>In the next examples we demonstrate how to produce predictions using multi-scale and/or multi-source data. By multi-scale data we imply covariates used for geostatistical mapping that are available at 2 or more (distinctly different) resolutions, but then cover the same area of interest (see also: RasterStack in the <a href="http://r-forge.r-project.org/projects/raster/" target="_blank">raster package</a>). Or in other words, we fit models and generate predictions with covariates at different resolutions, but with the same extent. In the case of the multi-source data, covariates can be of any scale, they can have a variable extent, and variable accuracy. As a general strategy, for multi-scale data we propose fitting a single model to combined covariates, while for the multi-source data we propose using data assimilation methods (merging of predictions). </p>
<table width="500" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Multi-scale versus multi-source approaches to  soil mapping.
  </caption>
  <tr>
    <th scope="col"><img src="Fig_multiscale_vs_multisource.png" alt="Fig_multiscale_vs_multisource.png" width="600" /></th>
  </tr>
</table>
<h3><a name="RK_multiscale" id="RK_multiscale"></a>Regression-kriging using multiscale data</h3>
<p><a href="http://plotkml.r-forge.r-project.org/eberg.html">Eberg&ouml;tzen</a> case study contains covariates at 100 and 25 m resolution. Both list of covariates have the same extent, which allows us to fit a model by still using the same <a href="fit.gstatModel-method.html">fit.gstatModel</a> function used above:</p>
<p class="R_code">&gt; eberg_grids &lt;- list(eberg_grid, eberg_grid25)<br />
&gt; unlist(sapply(eberg_grids, names))</p>
<p class="R_env">  [1] &quot;PRMGEO6&quot; &quot;DEMSRT6&quot; &quot;TWISRT6&quot; &quot;TIRAST6&quot; &quot;LNCCOR6&quot; &quot;DEMTOPx&quot; &quot;HBTSOLx&quot; &quot;TWITOPx&quot; &quot;NVILANx&quot;</p>
<p class="R_code">&gt; glm.formulaString3 = observedValue ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6+LNCCOR6+TWITOPx+NVILANx+ns(altitude, df=4)<br />
&gt; SNDMHT.m3 &lt;- fit.gstatModel(observations=eberg.geo, glm.formulaString3, covariates=eberg_grids, methodid=&quot;SNDMHT.t&quot;)</p>
<p><a href="fit.gstatModel-method.html">fit.gstatModel</a> function recognizes that the input to the model is a list of grids, and hence overlays points and grids in a loop. Note however, that the model fitting might be time consuming as the model needs to overlay points and grids, then merge resulting outputs in a loop.</p>
<p>The model above (<span class="R_code">SNDMHT.m3</span>) would be more difficult to use for prediction because it can happen that some land cover classes or soil mapping units do not contain any points. The model will thus not be able to predict values at new locations and R will report an error e.g. &quot;error in model... factor 'LNCCOR6' has new level(s)&quot;. Thus, a more robust way to generate predictions with multiscale covariates is to, again, convert the original predictors to principal components that are all continuous and independent. We first need to <a href="make.3Dgrid-method.html">downscale</a> all predictors from 100 m to 25 m resolution using the gdalwarp function:</p>
<p class="R_code">&gt; eberg_grid25@data &lt;- cbind(eberg_grid25@data, gdalwarp(eberg_grid, pixsize=eberg_grid25@grid@cellsize[1], GridTopology=eberg_grid25@bbox, resampling_method=&quot;cubicspline&quot;)@data)</p>
<p class="R_env">  Resampling 5 layers to &quot;+init=epsg:31467&quot; with grid cell size of 25 ...</p>
<p>So that we can derive soil predictive component via the <a href="spc-method.html">spc</a> function: </p>
<p class="R_code">&gt; formulaString2 &lt;- ~ TWITOPx+NVILANx+PRMGEO6+DEMSRT6+TWISRT6+TIRAST6<br />
  &gt; eberg_spc25 &lt;- spc(eberg_grid25, formulaString2)</p>
<table width="500" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: First four components derived using the combined 100 m and 25 m resolution data.
    Compare with the PCs derived using the 100 m resolution covariates only.
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_SPCs2_4.png" alt="Fig_eberg_SPCs2_4.png" width="500" /></th>
  </tr>
</table>
<p>which now allows us to model spatial distribution of sand content as a function of 14 components, at four times better resolution: </p>
<p class="R_code">&gt; glm.formulaString3 = as.formula(paste(&quot;observedValue ~ &quot;, paste(names(eberg_spc25@predicted), collapse=&quot;+&quot;), &quot;+ ns(altitude, df=4)&quot;))<br />
&gt; glm.formulaString3</p>
<p>  <span class="R_env">observedValue ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + <br />
PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + ns(altitude, df = 4)</span></p>
<p class="R_code">&gt; SNDMHT.m3 &lt;- fit.gstatModel(observations=eberg.geo, glm.formulaString3, covariates=eberg_spc25@predicted, methodid=&quot;SNDMHT.t&quot;) </p>
<p>Note that the resulting model (<span class="R_code">glm.formulaString3</span>) is only slightly better than the model we fitted using the 100 m resolution data. For mapping sand content in this study area, there seems to be not much gain in using 25 m resolution covariates considering the model significance. The difference is, however, that the predictions with 25 m resolution covariates  show four times better detail than the model with 100 m resolution grids (see small topographic features in the figure below), hence on the end there is a significant gain in using higher resolution data.</p>
<table width="500" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Comparison of predictions created using 25 m and 100 m resolution covariates. The missing values indicate areas of where the predictions were signficant (and hence were masked out).
  Values from 4% to 96% sand (blue to red).
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_comparison_25_100_m.png" alt="Fig_eberg_comparison_25_100_m.png" width="600" /></th>
  </tr>
</table>
<p>It is worth mentioning that preparation of the prediction locations for multiscale data is, however, slightly more complicated as the covariates had to be resampled from 100 m resolution to 25 m resolution. This operation can be time consuming and can lead to artifacts.  GSIF, by default, uses <a href="http://fwtools.maptools.org/">FWTools</a> (<a href="http://www.gdal.org/gdalwarp.html">gdalwarp</a> function) for resampling and reprojection of grids. FWTools is highly suitable for processing large data, and has many possibilities considering the downscaling algorithms. For example, via FWTools you can pass the resampling method <a href="http://www.gdal.org/gdalwarp.html">&quot;bicubicspline&quot;</a> (Keys, 1981), which is  effective for downscaling relatively smooth grids such as DEM or DEM parameters. The problem is that R needs to read and write spatial grids between FWTools, so that processing grids using FWTools can still take significant amount of time / RAM when working with large grids.</p>
<h3> <a name="RK_multisource" id="RK_multisource"></a>Merging multi-source data</h3>
<p>Imagine if the covariates for the <a href="http://plotkml.r-forge.r-project.org/eberg.html">Eberg&ouml;tzen</a> case study did not have the same extent, so that by overlaying the points we would loose significant parts of the callibration data as the covariates would not match. To avoid this, we can instead produce multiple models for each extent, then use some type of the data assimilation method to merge the predictions produced by various models. A sensitive approach to merge multiple predictions is to use the per-pixel accuracy to assign the weights, so that more accurate predictions receive more weight (<a href="http://dx.doi.org/10.1016/0016-7061(92)90002-O" target="_blank">Heuvelink and Bierkens, 1992</a>).</p>
<p>In practice, merging multi-source predictions takes only few more steps than in the case of multiscale data. First, we fit models using the two scales independently:</p>
<p class="R_code">&gt; formulaString.l &lt;- list(~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6, ~ DEMTOPx+TWITOPx+NVILANx)<br />
&gt; eberg_grids_spc &lt;- spc(eberg_grids, formulaString.l)</p>
<p class="R_env">  Converting PRMGEO6 to indicators...<br />
  Converting covariates to principal components...<br />
  Converting covariates to principal components...</p>
<p class="R_code">&gt; glm.formulaString.l &lt;- lapply(eberg_grids_spc, FUN=function(x){as.formula(paste(&quot;observedValue ~ &quot;, paste(names(x@predicted), collapse=&quot;+&quot;), &quot;+ ns(altitude, df=4)&quot;))})<br />
&gt; glm.formulaString.l</p>
<p class="R_env">  [[1]]<br />
  observedValue ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + <br />
  PC9 + PC10 + PC11 + ns(altitude, df = 4)<br />
&lt;environment: 0x000000001f65fff0&gt;</p>
<p class="R_env">[[2]]<br />
  observedValue ~ PC1 + PC2 + PC3 + ns(altitude, df = 4)<br />
&lt;environment: 0x0000000016d08830&gt;</p>
<p class="R_code">&gt; SNDMHT.ml &lt;- fit.gstatModel(observations=eberg.geo, glm.formulaString.l, lapply(eberg_grids_spc, slot, &quot;predicted&quot;), methodid=&quot;SNDMHT.t&quot;, rvgm=NULL)</p>
<p>Second, we can produce predictions using a list of covariates at different scales:</p>
<p class="R_code">&gt; new3D.ml &lt;- sapply(eberg_grids_spc, FUN=function(x){sp3D(x@predicted, stdepths=-0.025)})<br />
&gt; sd.ml &lt;- predict(SNDMHT.ml, predictionLocations=new3D.ml, nfold=2, verbose=TRUE, mask.extra=FALSE)</p>
<p class="R_env">Generating predictions using the trend model (KED method)...<br />
  [using universal kriging]</p>
<p class="R_env">Running 2-fold cross validation...<br />
|==============================================================================================| 100%</p>
<p class="R_env">Generating predictions using the trend model (KED method)...<br />
  [using universal kriging]</p>
<p class="R_env">Running 2-fold cross validation...<br />
|==============================================================================================| 100%</p>
<p>In the third step, we can sum up the predictions using weighted averaging. <a href="merge.html">merge</a> function will by default use results of cross-validation to scale the prediction variances, which are then used as weigths:</p>
<p class="R_code">&gt; sd.mlc &lt;- merge(sd.ml[[1]], sd.ml[[2]], silent=FALSE)</p>
<p class="R_env">Resampling 2 layers to &quot; +init=epsg:31467 +p &quot; with grid cell size: 25 ...<br />
  |==============================================================================================| 100%<br />
  Cross-validation RMSE (type = link):<br />
  observedValue_1 observedValue_2 <br />
0.5713          0.6340</p>
<p>In this case the predictions using the 100 m resolution data will receive a slightly higher weights. </p>
<table width="500" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Predictions using multi-source data (produced using the function merge).
Values from 4% to 96% sand (blue to red). 
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_merge_25_100_m.png" alt="Fig_eberg_merge_25_100_m.png" width="500" /></th>
  </tr>
</table>
<p>The results presented in figure above confirm that this method is equally valid and leads to almost identical predictions as when all covariates are used together (compare with the previous section). The final map again reflects mainly the patterns of the coarser predictors as these are more dominant. Assuming that the two predictions use completely independent covariates, the output should indeed be very similar to using a method that combines the two models. In fact, also in the case of regression-kriging, predictions using covariates are summed up to predictions using kriging (of residuals). In practice, the covariates might be correlated, and the models could overlap, so that the output of <a href="merge.html">merge</a> function might be biased. In fact, merging multi-source data is highly sensitive to type of models used to generate predictions and associated assumptions (for a discussion see e.g. <a href="http://web.gps.caltech.edu/~drf/misc/airs/maup_summary.pdf">Gotway and Young, 2002</a>; <a href="http://dx.doi.org/10.1198/106186007X179257" target="_blank">Gotway and Young, 2007</a>). Nevertheless, this approach is attractive for global soil mapping applications where predictions can be produced for areas of various extent, or using completely independent methods that can not be combined mathematically. For example, by following this approach, one could build models and produce predictions of soil properties for the whole world but at coarse resolution of e.g. 1 km, then use these predictions to fill-in the gaps for a regional model at resolution of 250 m or 100 m. It is important, however, to be aware that the covariates,  even if they are at different scales, are cross-correlated and hence this approach might lead to biased predictions (over or under estimation). </p>
<p>Finally, the propagated prediction uncertainty for merged predictions is more difficult to estimate as one can not simply take the mean of the prediction (kriging) variances. Instead, geostatistical simulations can be used to produce multiple realizations (see figure above). The propagated prediction error can then be derived by aggregating multiple realizations per pixel (see e.g. <a href="http://cran.r-project.org/web/packages/raster/" target="_blank">raster</a> package).</p>
<table width="100%" border="0" cellspacing="0" cellpadding="10">
  <tr>
    <th scope="col"><div align="left">
      <h2><a name="references" id="references"></a>References</h2>
    </div></th>
    <th scope="col">&nbsp;</th>
    <th scope="col"><div align="right"><a href="#top">^to top</a> </div></th>
  </tr>
</table>
<ol>
  <li>Bivand, R.S., Pebesma, E.J., and Gómez-Rubio, V., (2008) <a href="http://www.asdar-book.org/"><strong>Applied Spatial Data Analysis with R</strong></a>. Springer, 378 p.</li>
  <li>Böhner, J., McCloy, K. R. and Strobl, J. (Eds), (2006) <a href="http://saga-gis.org/en/about/references.html" target="_blank"><strong>SAGA &#8212; Analysis and Modelling Applications</strong></a>. Göttinger Geographische Abhandlungen, Heft 115. Verlag Erich Goltze GmbH, Göttingen, 117 pp.</li>
  <li> Gotway, C.A.,  Young, L.J., (2002) <a href="http://web.gps.caltech.edu/~drf/misc/airs/maup_summary.pdf" target="_blank">Combining Incompatible Spatial Data</a>. Journal of the American Statistical Association, <strong>97</strong> (458): 632-648. </li>
  <li>Gotway, C.A.,  Young, L.J., (2007) <a href="http://dx.doi.org/10.1198/106186007X179257" target="_blank">A Geostatistical Approach to Linking Geographically Aggregated Data From Different Sources</a>. Journal of Computational and Graphical Statistics, <strong>16</strong> (1), 115-135</li>
  <li> Hengl, T., Nikolic, M., MacMillan, R.A., (2012) <a href="http://dx.doi.org/10.1016/j.jag.2012.02.005" rel="nofollow">Mapping efficiency and information content</a>. International Journal of Applied Earth Observation and Geoinformation, special issue <em>Spatial Statistics Conference</em>, in press. </li>
  <li>Hijmans, R.J, Elith, J., (2012) <a href="http://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf">Species distribution modeling with R</a>. CRAN, Vignette for the dismo package, 72 p.</li>
  <li> Heuvelink, G.B.M., Bierkens, M.F.P. (1992) <a href="http://dx.doi.org/10.1016/0016-7061(92)90002-O" target="_blank">Combining soil maps with interpolations from point observations to predict quantitative soil properties</a>. Geoderma <strong>55</strong> (1-2): 1-15. </li>
  <li>Keys, R. (1981) <a href="http://dx.doi.org/10.1109/TASSP.1981.1163711">Cubic convolution interpolation for digital image processing</a>. <em>IEEE Transactions on Signal Processing, Acoustics, Speech, and Signal Processing</em> <strong>29</strong> (6): 1153–1160.</li>
  <li>Pebesma, E., (2003) <a href="http://www.gstat.org/gstat.pdf" target="_blank">Gstat user's manual</a>. Dept. of Physical Geography, Utrecht University, p. 100, www.gstat.org</li>
  <li>Skaggs, T. H., Arya, L. M., Shouse, P. J., Mohanty, B. P., (2001) <a href="https://www.soils.org/publications/sssaj/articles/65/4/1038" target="_blank">Estimating Particle-Size Distribution from Limited Soil Texture Data</a>. Soil Science Society of America Journal <strong>65</strong> (4): 1038-1044.</li>
  <li>Venables, W. N. and Ripley, B. D. (2002) <a href="http://www.stats.ox.ac.uk/pub/MASS4/"><strong>Modern Applied
Statistics with S</strong></a>. 4th edition.  Springer.</li>
</ol>
<p>&nbsp; </p>
<hr />
<p><a href="http://www.isric.org" target="_blank"><img src="ISRIC_logo.jpg" alt="ISRIC - World Soil Information" width="363" height="98" border="0" longdesc="http://www.isric.org" /></a></p>
</body>
</html>
