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
<p class="style1">Contact: <a href="http://www.wewur.wur.nl/popups/vcard.aspx?id=HENGL001" target="_blank">Tomislav Hengl</a><br />
  Last update:
  <!-- #BeginDate format:Am1 -->July 23, 2012<!-- #EndDate -->
</p>
<p>The purpose of this tutorial is to demonstrate major processing steps used within the GSIF framework to generate soil property and soil class maps from point data and with the help of multi-scale covariates. The GSIF (R package) <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. To learn more about the Global Soil Information Facilities (GSIF), visit the <a href="http://www.isric.org/projects/global-soil-information-facilities-gsif" target="_blank">main project page</a>. See also the complete list of <strong><a href="00Index.html">functions</a></strong> available via the GSIF package.</p>
<p>Download the tutorial as <a href="tutorial_eberg.R">R script</a>. </p>
<hr />
<h2>Loading the data and  data screening</h2>
<p>For demonstration purposes we use the <a href="http://plotkml.r-forge.r-project.org/eberg.html">Eberg&ouml;tzen</a> case study, which   has been used by the SAGA GIS development team. To start the exercise, first install and load all required packages:</p>
<p class="R_code">&gt; library(plotKML)<br />
&gt; library(GSIF)</p>
<p class="R_env">Loading required package: RCurl<br />
Loading required package: bitops<br />
GSIF version 0.2-1 (2012-07-19)<br />
URL: http://gsif.r-forge.r-project.org/</p>
<p>then load the input data:</p>
<p class="R_code">&gt; data(eberg)<br />
  &gt; data(eberg_grid)<br />
  &gt; data(eberg_grid25)</p>
<p>The <a href="http://plotkml.r-forge.r-project.org/eberg.html">Eberg&ouml;tzen</a> case study is a 10 by 10 km large case study in central Germany. The eberg object contains soil profile observations, and <span class="R_code">eberg_grid</span> and <span class="R_code">eberg_grid25</span> contains covariate layers at 100 and 25 m resolution. Before we proceed with running geostatistical analysis, we can run some initial data screening. First, we look at the data structure: </p>
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
<p>This shows that soil texture fractions <span class="R_code">SNDMHT</span>, <span class="R_code">SLTMHT</span> and <span class="R_code">CLYMHT</span> have been sampled at fixed depths (10-30, 30-50, 50-70, 70-90), and <span class="R_code">soiltype</span> comprises 13 soil types. Let us see some summary info for sand content:</p>
<p class="R_code">&gt; summary(eberg$SNDMHT_A)</p>
<p class="R_env">Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's <br />
  4.062  18.800  20.000  30.260  42.300  92.500   3.000 </p>
<p class="R_code">&gt; library(StatDA)<br />
  &gt; par(mar=c(2.5,2.5,0.5,0.5), oma=c(0,0,0,0))<br />
  &gt; edaplot(eberg$SNDMHT_A[!is.na(eberg$SNDMHT_A)], H.freq=TRUE, box=FALSE, S.pch=3, S.cex=0.5, D.lwd=1.5, P.ylab=&quot;&quot;, P.log=FALSE, P.logfine=c(5,10), P.main=&quot;&quot;, P.xlab=&quot;&quot;, B.pch=3, B.cex=0.5)</p>
<table width="650" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Histogram for sand content.
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_hist_SNDMHT.png" alt="Fig_eberg_hist_SNDMHT.png" width="634" height="638" /></th>
  </tr>
</table>
<p>From the plot above, we can also observe that many numbers are in fact overlapping. It seems that there is only a clusters of<br />
  values possible for <span class="R_code">SNDMHT</span>. Going 
  back to the origin of this dataset, we can notice that the sand, silt and clay values have 
  been determined by using the so-called texture by hand method, i.e. via texture classes. The literature (<a href="https://www.soils.org/publications/sssaj/articles/65/4/1038" target="_blank">Skaggs
  et al., 2001</a>) reports that the this technique can be used to determine the content of soil
earth fractions only to an accuracy of ±5–10%. This means that we should not plan to map any of the texture fractions to a precision better than ±5% (detection limit), because it would exceed the measurement error.</p>
<p>These is a relatively large set (<span class="R_code">3670</span> points), so it might be a good idea to subset it (e.g. take only 30% of samples) to speed up processing:</p>
<p class="R_code">&gt; eberg.xy &lt;- eberg[runif(nrow(eberg)) &lt; .3,]<br />
  &gt; coordinates(eberg.xy) &lt;- ~X+Y<br />
  &gt; proj4string(eberg.xy) &lt;- CRS(&quot;+init=epsg:31467&quot;)<br />
  &gt; gridded(eberg_grid) &lt;- ~x+y<br />
  &gt; proj4string(eberg_grid) &lt;- CRS(&quot;+init=epsg:31467&quot;)<br />
  &gt; gridded(eberg_grid25) &lt;- ~x+y<br />
  &gt; proj4string(eberg_grid25) &lt;- CRS(&quot;+init=epsg:31467&quot;)</p>
<p>To assess how representative and consistent 
  is soil profile data, we can run some basic analysis of the point geometry and<br />
  then overlap the points with predictors to see how well are the environmental features 
  represented. First, we can test if the points represent the geographical space:</p>
<p class="R_code">&gt; libary(spatstat) <br />
&gt; eberg.ppp &lt;- ppp(x=coordinates(eberg.xy)[,1], y=coordinates(eberg.xy)[,2], marks=eberg.xy$zinc, window=mg_owin)<br />
&gt; summary(nndist(eberg.ppp))</p>
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
<table width="650" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: The results of the Complete Spatial Randomness test (based on spatstat) for the Eberg&ouml;tzen case study.
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_CRS_test.png" alt="Fig_eberg_CRS_test.png" width="678" height="657" /></th>
  </tr>
</table>
<p>which shows that the point samples do not exactly satisfy the Complete Spatial Randomness test &#8212; samples at larger distances have a lower spatial density than a completely random design, which usually means that large parts of the case study are undersampled (see also figure below).</p>
<p>To see how representative is the point data considering the coverage of feature space, we use the maxent function, available via the dismo package:</p>
<p class="R_code">&gt; me.eberg &lt;- MaxEnt(occurrences=eberg.ppp, covariates=eberg_grid)<br />
  &gt; par(mfrow=c(1,2), mar=c(0.5,0.5,0.5,0.5), oma=c(0,0,0,0))<br />
  &gt; image(as(me.eberg@predicted, &quot;SpatialPixelsDataFrame&quot;), col=rev(heat.colors(25)), xlab=&quot;&quot;, ylab=&quot;&quot;)<br />
  &gt; points(me.eberg@occurrences, pch=&quot;+&quot;, cex=.7)<br />
  &gt; image(me.eberg@sp.domain, col=&quot;grey&quot;, xlab=&quot;&quot;, ylab=&quot;&quot;)</p>
<table width="800" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Niche analysis (sampling likelihood in feature space) based on MaxEnt (left; dark red indicates higher values) and areas fully represented by the current samples following the cross-validation (right).
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_MaxEnt_test.png" alt="Fig_eberg_MaxEnt_test.png" width="857" height="533" /></th>
  </tr>
</table>
<p>Which shows that, some parts of the study area (higher elevations, specific land cover types) have been systematically omitted from sampling. The map on the right shows which pixels (grey) are actually valid to run predictions as the probability of occurrence of sampling points at these locations is statistically significant. This is nothing that should worry us too much, but something we need to be aware when doing the interpretation of produced maps. To learn  more about MaxEnt, refer to the <a href="http://cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf" target="_blank">dismo package vignette</a>. </p>
<hr />
<h2>Model building and spatial predictions </h2>
<p>Any geostatistical mapping boils down to two major processes: (<em>a</em>) model fitting (or model estimation), and (<em>b</em>) generation of spatial predictions. Also in GSIF package geostatistical mapping is implemented via two main functions: <a href="fit.gstatModel-method.html">fit.gstatModel</a> and predict. In the following sections we demonstrate how to fit models and generate predictions using soil profile data and a diversity of covariate layers, typical for soil mapping applications. </p>
<p>Before we build a regression-kriging model to generate 3D predictions, however, we need to prepare the soil observations to a format suitable for fitting of 3D models. We start by converting the <span class="R_code">eberg</span> data frame to class SoilProfileCollection-class (used by the <a href="http://cran.r-project.org/package=aqp">aqp</a> package) and the <a href="http://gsif.r-forge.r-project.org/geosamples-class.html">geosamples-class</a> (used by the GSIF package):</p>
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
<p>GSIF package by default works with geosamples as the main class for point observations. Convertion of the SPC class data to geosamples is straight forward:</p>
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
..@ data    :'data.frame':    80740 obs. of  13 variables:<br />
.. ..$ observationid   : chr [1:80740] NA NA NA NA ...<br />
.. ..$ sampleid        : chr [1:80740] &quot;id0001&quot; &quot;id0002&quot; &quot;id0003&quot; &quot;id0004&quot; ...<br />
.. ..$ longitude       : num [1:80740] 10.2 10.2 10.2 10.2 10 ...<br />
.. ..$ latitude        : num [1:80740] 51.6 51.6 51.6 51.6 51.5 ...<br />
.. ..$ locationError   : num [1:80740] NA NA NA NA NA NA NA NA NA NA ...<br />
.. ..$ TimeSpan.begin  : POSIXct[1:80740], format:  ...<br />
.. ..$ TimeSpan.end    : POSIXct[1:80740], format:  ...<br />
.. ..$ altitude        : num [1:80740] 0 0 0 0 0 0 0 0 0 0 ...<br />
.. ..$ altitudeMode    : chr [1:80740] &quot;relativeToGround&quot; &quot;relativeToGround&quot; &quot;relativeToGround&quot; &quot;relativeToGround&quot; ...<br />
.. ..$ volume          : num [1:80740] 2 2 2 2 2 2 2 2 2 2 ...<br />
.. ..$ observedValue   : chr [1:80740] &quot;G&quot; &quot;L&quot; &quot;L&quot; &quot;G&quot; ...<br />
.. ..$ methodid        : Factor w/ 6 levels &quot;CLYMHT&quot;,&quot;SLTMHT&quot;,..: 5 5 5 5 5 5 5 5 5 5 ...<br />
.. ..$ measurementError: num [1:80740] NA NA NA NA NA NA NA NA NA NA ...</p>
<p class="R_code">&gt; levels(eberg.geo@data$methodid)</p>
<p class="R_env">  [1] &quot;CLYMHT&quot;   &quot;SLTMHT&quot;   &quot;SNDMHT&quot;   &quot;SNDMHT.t&quot; &quot;soiltype&quot;<br />
  [6] &quot;TAXGRSC&quot;</p>
<p>Geosamples-class can be considered the most plain (standard) format for any space-time observations, and main advantage of using this class in R is that it can be easily manipulated and converted to spatial and aqp classes. The column names in the <span class="R_code">@data</span> slot correspond to the tag names used in the KML schema, which makes it easier to export such data to some GIS or Google Earth. The dissadvantage of using geosamples is that the object size can easily grow because each observation gets a separate row.</p>
<p>Prior to geostatistical modelling, it is also probably a good idea to convert all covariates to independent components. This way, it will be easier to subset to the optimal number of predictors during the regression analysis. PCA helps also reducing the prediction bias, assuming that the covariates are cross-correlated. A wrapped function <a href="spc-method.html">spc</a> will convert all factor variables to indicators and run PCA on a stack of grids: </p>
<p class="R_code">&gt; formulaString &lt;- ~ PRMGEO6+DEMSRT6+TWISRT6+TIRAST6<br />
&gt; eberg_spc &lt;- spc(eberg_grid, formulaString)</p>
<p class="R_env">Converting PRMGEO6 to indicators...<br />
  Converting covariates to principal components...</p>
<p>The output is a stack of independent components, all numeric and all scaled around 0 value. To see which inputs define any of the components, we can look at the rotation table or plot the images: </p>
<p class="R_code">&gt; eberg_spc@pca$rotation<br />
&gt; pal = rev(rainbow(65)[1:48])<br />
  &gt; rd = range(eberg_spc@predicted@data[,1], na.rm=TRUE)<br />
&gt; spplot(eberg_spc@predicted[1:4], at=seq(rd[1], rd[2], length.out=48), col.regions=pal)</p>
<table width="600" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: First four components derived using the eberg_grid data.
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_SPCs1_4.png" alt="Fig_eberg_SPCs1_4.png" width="587" height="537" /></th>
  </tr>
</table>
<h3>1. Soil properties </h3>
<p>In the GSIF package, regression-kriging model can be fitted at once by using the <a href="fit.gstatModel-method.html">fit.gstatModel</a> function. First, we need to specify the soil property ~ covariates model: </p>
<p class="R_code">&gt; glm.formulaString = as.formula(paste(&quot;observedValue ~ &quot;, paste(names(eberg_spc@predicted), collapse=&quot;+&quot;), &quot;+ ns(altitude, df=4)&quot;))<br />
&gt; glm.formulaString</p>
<p class="R_env">  observedValue ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + <br />
PC9 + PC10 + PC11 + ns(altitude, df = 4)</p>
<p>In other words the observed values will be modelled as a function of PCs and altitude (natural splines). 3D GLM-kriging model can now be fitted at once: </p>
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
<p>These result show that the model is significant, both the GLM and the variogram. Note however that it we do not actually have values of PCs at different depths (in fact, most of PCs relate only to the surface), so that many values of covariates are basically copied to all depths. This does not represent any problem for the GLM modelling, however, you should be aware that, because values of covariates are fixed with the different depths, the 3D patterns will be mainly controlled by the surface patterns.</p>
<p>Now that we have fitted a <a href="gstatModel-class.html">gstatModel</a>, we can basically generate predictions and estimate the associated uncertainty at any depth. In the last step, we need to prepare the 3D prediction locations i.e. grid cells that need to be mapped. In GSIF package, this can be done by using the <a href="make.3Dgrid-method.html">sp3D</a> function:</p>
<p class="R_code">&gt; new3D &lt;- sp3D(eberg_spc@predicted)</p>
<p>This will prepare the existing SPCs for 3D predictions (6 standard depths) by adding the altitude column:</p>
<p class="R_code">&gt; str(new3D[[1]]@grid)</p>
<p class="R_env">  Formal class 'GridTopology' [package &quot;sp&quot;] with 3 slots<br />
  ..@ cellcentre.offset: Named num [1:3] 3.57e+06 5.71e+06 -2.50e-02<br />
  .. ..- attr(*, &quot;names&quot;)= chr [1:3] &quot;longitude&quot; &quot;latitude&quot; &quot;altitude&quot;<br />
  ..@ cellsize         : Named num [1:3] 1e+02 1e+02 5e-02<br />
  .. ..- attr(*, &quot;names&quot;)= chr [1:3] &quot;longitude&quot; &quot;latitude&quot; &quot;altitude&quot;<br />
  ..@ cells.dim        : Named int [1:3] 100 100 1<br />
  .. ..- attr(*, &quot;names&quot;)= chr [1:3] &quot;longitude&quot; &quot;latitude&quot; &quot;altitude&quot; </p>
<p>So that we can now run prediction for each standard depth in a loop: </p>
<p class="R_code">&gt; sd.l &lt;- lapply(new3D, FUN=function(x){predict(SNDMHT.m, predictionLocations=x, nfold=0)})</p>
<p class="R_env">Generating predictions using the trend model (KED method)...<br />
  [using universal kriging]</p>
<p>This operation can take time depending on the size of the grids and number of 3D points used to generate predictions. Finally, we can prepare the produced predictions and export them as <a href="GlobalSoilMap-class.html">GlobalSoilMap-class</a> object. First, we back-transform the predictions to the 0-100% scale:</p>
<p class="R_code">&gt; for(j in 1:length(sd.l)){ sd.l[[j]]@predicted$observedValue &lt;- exp(sd.l[[j]]@predicted$observedValue)/(1+exp(sd.l[[j]]@predicted$observedValue))*100 }</p>
<p>Next, we can reproject the produced predictions to the WGS84 geographical coordinates: </p>
<p class="R_code">&gt; p = get(&quot;cellsize&quot;, envir = GSIF.opts)[2]<br />
  &gt; s = get(&quot;stdepths&quot;, envir = GSIF.opts)<br />
  &gt; s</p>
<p class="R_env">  [1] -0.025 -0.075 -0.225 -0.450 -0.800 -1.500</p>
<p class="R_code">&gt; sd.ll &lt;- sapply(1:length(sd.l), FUN=function(x){make.3Dgrid(sd.l[[x]]@predicted[3:4], pixelsize=p, stdepths=s[x])})</p>
<p class="R_env">  Resampling 2 layers to +proj=longlat +datum=WGS84 with grid cell size of 0.000833333333333333 ...</p>
<p>The object <span class="R_code">sd.ll</span> can now be saved as a GlobalSoilMap object:</p>
<p class="R_code">&gt; SNDMHT.gsm &lt;- GlobalSoilMap(varname=&quot;SNDMHT&quot;, sd.ll, period=c(&quot;1999-02-01&quot;, &quot;2001-07-01&quot;))<br />
<span class="R_code">&gt; save(SNDMHT.gsm, file=&quot;SNDMHT.rda&quot;, compress=&quot;xz&quot;)</span>  </p>
<p>and visualized  in Google Earth by using the plotKML package.</p>
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
<table width="600" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Predicted sand content using 3D GLM-kriging as visualized in Google Earth (<a href="eberg_SNDMHT_6depths.kmz">kmz</a>). 
  </caption>
  <tr>
    <th scope="col"><a href="eberg_SNDMHT_6depths.kmz"><img src="Fig_GSIF_predicted_layers.jpg" alt="Fig_GSIF_predicted_layers.jpg" width="600" border="0" /></a></th>
  </tr>
</table>
<h3>2. Soil Classes </h3>
<p>GSIF package also provides functionality for pedometric mapping of soil classes. Soil types can be mapped using a wrapper function   <a href="spfkm.html">spfkm</a>. This will   run supervised fuzzy <em>k</em>-means using a list of covariates layers provided as <span class="R_code">&quot;SpatialPixelsDataFrame&quot;</span> object, and optional classe centres and class variances. As in the case of continuous/numeric variables, the process consists of model fitting and predictions, which is in this case wrapped into a single function: </p>
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
<p>The output is an object of class <a href="SpatialMemberships-class.html">SpatialMemberships</a>. In the case above,  the class centres and variances were not specified, hence spfkm tries to estimate them via the multinomial logistic regression model (nnet::multinom). They could have also been passed via the <span class="R_code">class.c</span> and <span class="R_code">class.sd</span> arguments i.e. based on the expert judgement. To see the actual class centres use:</p>
<p class="R_code">&gt; eberg_sm@class.c[1,]</p>
<p class="R_env"> PC1        PC2        PC3        PC4        PC5        PC6        PC7        PC8 <br />
3.1607835  1.7982977  0.2260035 -0.5133303 -0.9440280 -0.2689071  0.4962124  1.9947157 <br />
PC9       PC10 <br />
0.8995294 -0.3908947</p>
<p class="R_code">&gt; row.names(eberg_sm@class.c) </p>
<p class="R_env">[1] &quot;A&quot;  &quot;B&quot;  &quot;D&quot;  &quot;G&quot;  &quot;Hw&quot; &quot;L&quot;  &quot;N&quot;  &quot;Q&quot;  &quot;R&quot;  &quot;S&quot;  &quot;Z&quot;</p>
<table width="600" border="0" cellspacing="2" cellpadding="4">
  <caption class="caption" align="bottom">
    Fig: Predicted soil types for the Eberg&ouml;tzen case study.
  See <a href="spfkm.html">spfkm</a> for more detatils.
  </caption>
  <tr>
    <th scope="col"><img src="Fig_eberg_Soiltypes_spfkm.png" alt="Fig_eberg_Soiltypes_spfkm.png" width="514" height="476" /></th>
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
<p>Note that intercept needs to be taken out, so that the best predictor of the sand content for some soil type is basically the mean value of the sand for that soil type (for example class B is expected to have an average sand content of about 50.1%). If you compare the model based on soil classes and model fitted in the previous section (<span class="R_code">SNDMHT.m</span>), you can see that fitting data using a 3D model results in a slightly better fit. Nevertheless, soil classes are in this case study significant estimators of sand content. </p>
<hr />
<h2>Multiscale modelling</h2>
<p>In the next examples we demonstrate how to produce predictions using multiscale / multisource data. Multiscale data implies that covariates are available at 2 or more distinctly different resolutions, but then cover the same area of interest (different resolutions, same extent). Multisource data is a more general case of covariate data. Multisource data can have: (1) variable scale, (2) variable extent, (3) variable accuracy. </p>
<h3>1. All covariates available at all grid nodes</h3>
<p>&nbsp;</p>
<h3>2. Merging multisource data </h3>
<p>&nbsp;</p>
<hr />
<p><a href="http://www.isric.org" target="_blank"><img src="ISRIC_logo.jpg" alt="ISRIC - World Soil Information" width="363" height="98" border="0" longdesc="http://www.isric.org" /></a></p>
</body>
</html>
