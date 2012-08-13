
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
.style1 {font-size: x-small}
.R_code {font-family:"Courier New", Courier, monospace;
font-style:italic;
font-size: x-small;
}
-->
    </style>
</head>

<body>

<!-- R-Forge Logo -->
<table border="0" width="100%" cellspacing="0" cellpadding="0">
<tr><td>
<a href="http://r-forge.r-project.org/"><img src="<?php echo $themeroot; ?>imagesrf/logo.png" border="0" alt="R-Forge Logo" /> </a> </td> </tr>
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
  <p>This package contains tools and procedures to handle soil data and produce gridded soil property maps to support the global soil data initivatives such as the GlobalSoilMap.net project.</p>
</div>
<table border="0" cellspacing="0" cellpadding="10">
  <tr>
    <td><div align="center"><a href="http://globalsoilmap.net/content/1-degree-tiles-world"><img src="Fig_1degree_tiles_world.thumbnail.png" alt="One degree tiles (world map)" width="220" height="120" border="0" /></a></div></td>
    <td><a href="http://soilprofiles.org"><img src="opensoilprofiles.gif" alt="Open Soil Profiles" width="320" height="100" border="0" longdesc="http://soilprofiles.org" /></a></td>
  </tr>
</table>
<p class="style1">Contact: <a href="http://www.wewur.wur.nl/popups/vcard.aspx?id=HENGL001" target="_blank">Tomislav Hengl</a></p>
<p class="style1">Contributions by: Bas Kempen, Dylan Beaudette, Gerard B.M. Heuvelink, Reuter I. Hannes, Brendan Malone, Pierre Roudier ... </p>
<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>. See the complete list of <strong><a href="00Index.html">functions</a></strong> available in this package. GSIF provides access to several case studies. A tutorial to analyze soil property and soil class data with the Ebergotzen data is available <strong><a href="tutorial_eberg.php">here</a></strong>. </p>

<p><strong>Installation:</strong></p>
<p>To install this package from R-forge use (works only on<strong> &gt;= R 2.14!</strong>):</p>
<p class="R_code">&gt; install.packages(c(&quot;RCurl&quot;, &quot;XML&quot;, &quot;rgdal&quot;, &quot;raster&quot;, &quot;sp&quot;, &quot;aqp&quot;,  &quot;mda&quot;, &quot;gstat&quot;, &quot;plotKML&quot;, &quot;dismo&quot;))<br />
&gt; install.packages(&quot;GSIF&quot;, repos=c(&quot;http://R-Forge.R-project.org&quot;)) </p>
<p>Alternatively, you can install the most recent snapshot of the package directly from the source by using e.g.:</p>
<pre class="R_code">&gt; install.packages(c(&quot;RCurl&quot;, &quot;XML&quot;, &quot;rgdal&quot;, &quot;raster&quot;, &quot;sp&quot;, &quot;aqp&quot;,  &quot;mda&quot;, &quot;gstat&quot;, &quot;plotKML&quot;, &quot;dismo&quot;))
&gt; download.file(&quot;http://gsif.r-forge.r-project.org/GSIF_0.2-2.tar.gz&quot;, &quot;GSIF_0.2-2.tar.gz&quot;)
&gt; system(&quot;R CMD INSTALL GSIF_0.2-2.tar.gz&quot;) </pre>
<p>GSIF package extensively uses a number of external software, hence it is highly recommended that, prior to starting GSIF, you first download and install:</p>
<ol>
  <li><a href="http://fwtools.maptools.org" target="_blank">FWTools</a> (<em>required</em>) &#8212; this software is called by several functions (e.g. <a href="make.3Dgrid-method.html">make.3Dgrid</a>);</li>
  <li><a href="http://www.saga-gis.org" target="_blank">SAGA GIS</a> (<em>optional</em>) &#8212; this software is highly recommend but not required;</li>
  <li><a href="http://www.python.org/getit/" target="_blank">Python</a> (<em>optional</em>) &#8212; this allows processing of large data;</li>
</ol>
<p>Read more:  installation and first steps with <strong><a href="http://plotkml.r-forge.r-project.org/">plotKML</a></strong>.</p>
<p><strong>News:</strong></p>
<ul>
  <li>Aug 2012: added examples of how to <a href="http://plotkml.r-forge.r-project.org/tutorial.php">visualize</a> various (spatial) soil data; </li>
  <li>July 2012: added functionality for <a href="tutorial_eberg.php">3D regression-kriging and soil-class mapping</a>; </li>
  <li>Apr 2012: first meeting of the package development team at the <a href="http://www.pedometrics.org/dsm_oz/" target="_blank">DSM conference</a> in Sydney; </li>
  <li>Mar 2012: first version of the package on R-forge; </li>
</ul>
<hr />
<p><a href="http://www.isric.org" target="_blank"><img src="ISRIC_logo.jpg" alt="ISRIC - World Soil Information" width="363" height="98" border="0" longdesc="http://www.isric.org" /></a></p>
</body>
</html>
