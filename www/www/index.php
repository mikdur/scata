<?php
require_once("../includes/includes.php");

/* Hantera menyer
 * En visas om användare inloggad.
 * En annan om användare inte inloggad.
 *
 * Sidor som kräver att användare är inloggad har tre värden: Namn, del av URL och 1.
 * Sidor som inte kräver inloggning har samma värden, med undantag för 0 istället för 1.
 * Ska ett menyval växla mellan två sidor är arrayen fyra platser stor, med sidan för utloggad användare först.
 */
$menu_pages = Array();
$menu_pages[] = Array("Home", "home", "Home", "home");
$menu_pages[] = Array("Jobs", "jobs", 1);
$menu_pages[] = Array("Datasets", "datasets", 1);
$menu_pages[] = Array("Reference Sequences", "refseqs", 1);
$menu_pages[] = Array("Tag Sets", "tagsets", 1);
$menu_pages[] = Array("Parameter Sets", "parametersets", 1);
$menu_pages[] = Array("Files", "upload", 1);
$menu_pages[] = Array("Login", "login", "Logout", "logout");
$menu_pages[] = Array("Register", "register", 0);
$menu_pages[] = Array("Help", "help", "Help", "help");
$menu_pages[] = Array("Lost password", "lost_password", -1);

/**
 * Skriv ut meny
 */
function print_menu($menu_pages, $active)
{
	$user = new User();
	$page_num = 0;

	echo "<ul>";

	for ($i = 0; $i < count($menu_pages); $i++)
	{
		$page = "";
		$class = "inactive";
		if (count($menu_pages[$i]) == 3)
		{
			$page_num = 1;
			if ($user->is_logged_in())
			{
				if ($menu_pages[$i][2] == 1)
				{
					$page = "<a href=\"./?p={$menu_pages[$i][1]}\">{$menu_pages[$i][0]}</a>";
				}
				elseif ($menu_pages[$i][2] == 0)
				{
					$page = "{$menu_pages[$i][0]}</li>";
					$class = "grey";
				}
			}
			else
			{
				if ($menu_pages[$i][2] == 0)
				{
					$page = "<a href=\"./?p={$menu_pages[$i][1]}\">{$menu_pages[$i][0]}</a>";
				}
				elseif ($menu_pages[$i][2] == 1)
				{
					$page = "{$menu_pages[$i][0]}";
					$class = "grey";
				}
			}
		}
		elseif (count($menu_pages[$i]) == 4)
		{
			$page_num = ($user->is_logged_in() ? 1 : 0) * 2;
			$page = "<a href=\"./?p={$menu_pages[$i][($page_num+1)]}\">{$menu_pages[$i][$page_num]}</a>";
		}
		if ($menu_pages[$i][$page_num] == $active)
		{
			$class = "active";
		}
		if (strcmp($page,""))
		{
			echo "<li class=\"{$class}\">{$page}</li>\n";
		}
	}
	echo "</ul>";
}


/**
 * Kollar om användare begärt sida
 * Krävs inloggning verifieras det innan sidan visas
 */
function get_page_to_show($requested_page, $menu_pages)
{
	$user = new User();
	$page = "";

	if ($requested_page == "")
	{
		if ($user->is_logged_in())
		{
			return "home";
		}
		else
		{
			return "home";
		}
	}

	for ($i = 0; $i < count($menu_pages); $i++)
	{
		if (count($menu_pages[$i]) == 3)
		{
			if ($menu_pages[$i][1] == $requested_page)
			{
				if (($menu_pages[$i][2] == 0 && !$user->is_logged_in()) || ($menu_pages[$i][2] == 1 && $user->is_logged_in()) || ($menu_pages[$i][2] == -1))
				{
					$page = $menu_pages[$i][1];
					break;
				}
			}
		}
		elseif (count($menu_pages[$i]) == 4)
		{
			if ($menu_pages[$i][($user->is_logged_in() ? 1 : 0) * 2 + 1] == $requested_page)
			{
				$page = $menu_pages[$i][($user->is_logged_in() ? 1 : 0) * 2 + 1];
				break;
			}
		}
	}
	
	if ($page == "")
	{
		$page = "home";
	}
	
	return $page;
}

$page = get_page_to_show((isset($_REQUEST['p']) ? $_REQUEST['p'] : ""), $menu_pages);
?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.1//EN" "http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd"> 
<html xml:lang="sv" xmlns="http://www.w3.org/1999/xhtml"> 
<head> 
	<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
	<title>SCATA</title>
	<script type="text/javascript" src="./js/jquery-1.4.2.js"></script>
	<script type="text/javascript" src="./js/default.js"></script>
	<script type="text/javascript" src="./js/functions.js"></script>
<!-- production -->
<script type="text/javascript" src="./plupload/js/plupload.full.min.js"></script>


<!-- debug 
<script type="text/javascript" src="./js/moxie.js"></script>
<script type="text/javascript" src="./js/plupload.dev.js"></script>
-->

	<link rel="stylesheet" type="text/css" href="./css/style.css" />
</head>
<body>
<div id="container">
	<div id="navigation">
		<div id="nav-top-bar">
		</div>
		<div id="nav-logo">
			<img src="./images/menu/slu_1.gif" id="nav-logo-1" class="logo" />
			<img src="./images/menu/slu_2.gif" id="nav-logo-2" class="logo" />
		</div>
		<div id="nav-contents">
			<div id="nav-contents-filler">
				<h1>SCATA</h1>
				<span id="nav-contents-filler-text">
					Sequence Clustering and Analysis of Tagged Amplicons<br />
					Department of Forest Mycology and Pathology
				</span>
			</div>
			<div id="nav-scata-logo">
			     <img src="./images/menu/scata77.png" id="nav-scata-logo1" />
			</div>
			<div id="nav-contents-links">
				<?php
				print_menu($menu_pages, $page);
				?>
			</div>
		</div>
		<div id="nav-bottom-bar">
		</div>
	</div>
	<div id="contents">
		<?php
		require_once("../includes/pages/{$page}.php");
		?>
	</div>
</div>
</body>
</html>
