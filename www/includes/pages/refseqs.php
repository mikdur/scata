	<?php
	if (isset($_GET['do']))
	{
		switch ($_GET['do'])
		{
			case "add" :
				require_once("../includes/pages/refseqs/add.php");
				break;
			case "list" :
				require_once("../includes/pages/refseqs/list.php");
				break;
			case "delete" :
				require_once("../includes/pages/refseqs/delete.php");
				break;
			default :
				require_once("../includes/pages/refseqs/list.php");
				break;
		}
	}
	else
	{
		require_once("../includes/pages/refseqs/list.php");
	}
	?>