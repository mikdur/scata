	<?php
	if (isset($_GET['do']))
	{
		switch ($_GET['do'])
		{
			case "add" :
				require_once("../includes/pages/tagsets/add.php");
				break;
			case "list" :
				require_once("../includes/pages/tagsets/list.php");
				break;
			case "delete" :
				require_once("../includes/pages/tagsets/delete.php");
				break;
			default :
				break;
		}
	}
	else
	{
		require_once("../includes/pages/tagsets/list.php");
	}
	?>