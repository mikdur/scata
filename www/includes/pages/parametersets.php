<?php
if (isset($_GET['do']))
{
	switch ($_GET['do'])
	{
		case "add" :
			require_once("../includes/pages/parametersets/add.php");
			break;
		case "list" :
			require_once("../includes/pages/parametersets/list.php");
			break;
		case "details" :
			require_once("../includes/pages/parametersets/details.php");
			break;
		case "delete" :
			require_once("../includes/pages/parametersets/delete.php");
			break;
		default :
			require_once("../includes/pages/parametersets/list.php");
			break;
	}
}
else
{
	require_once("../includes/pages/parametersets/list.php");
}
?>