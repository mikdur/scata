	<?php
	if (isset($_GET['do']))
	{
		switch ($_GET['do'])
		{
			case "add" :
				require_once("../includes/pages/jobs/add.php");
				break;
			case "list" :
				require_once("../includes/pages/jobs/list.php");
				break;
			case "details" :
				require_once("../includes/pages/jobs/details.php");
				break;
			case "delete" :
				require_once("../includes/pages/jobs/delete.php");
				break;
			case "download" :
				require_once("../includes/pages/jobs/download.php");
				break;
			default :
				require_once("../includes/pages/jobs/list.php");
				break;
		}
	}
	else
	{
		require_once("../includes/pages/jobs/list.php");
	}
	?>