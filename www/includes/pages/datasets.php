	<?php
	if (isset($_GET['do']))
	{
		switch ($_GET['do'])
		{
			case "add" :
				require_once("../includes/pages/datasets/add.php");
				break;
			case "list" :
				require_once("../includes/pages/datasets/list.php");
				break;
			case "delete" :
				require_once("../includes/pages/datasets/delete.php");
				break;
			default :
				require_once("../includes/pages/datasets/list.php");
				break;
		}
	}
	else
	{
		require_once("../includes/pages/datasets/list.php");
	}
	?>