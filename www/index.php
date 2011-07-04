
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
.style2 {font-size: small}
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
  <p>Global Soil Information Facilities - tools and procedures to handle soil data and produce gridded soil property maps.</p>
</div>
<table border="0" cellspacing="0" cellpadding="10">
  <tr>
    <td><div align="center"><a href="http://globalsoilmap.org/content/1-degree-tiles-world"><img src="Fig_1degree_tiles_world.thumbnail.png" alt="One degree tiles (world map)" width="220" height="120" border="0" /></a></div></td>
    <td><img src="opensoilprofiles.gif" alt="Open Soil Profiles" width="320" height="100" longdesc="http://soilprofiles.org" /></td>
  </tr>
</table>
<p class="style1">Contact: <a href="http://www.wewur.wur.nl/popups/vcard.aspx?id=HENGL001" target="_blank">Tomislav Hengl</a></p>

<p> The <strong>project summary page</strong> you can find <a href="http://<?php echo $domain; ?>/projects/<?php echo $group_name; ?>/"><strong>here</strong></a>.</p>

<hr />
<p><img src="ISRIC_logo.jpg" alt="ISRIC - World Soil Information" width="297" height="72" longdesc="http://www.isric.org" /></p>
</body>
</html>
