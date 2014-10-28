<?php
/**
 * Flytta uppladdad fil
 * @param string Namn på uppladdad fil
 * @param string Nytt namn på fil
 * @param int Typ av fil
 * @return boolean True om det gick att flytta filen, false annars.
 */
function move_file($name, $new_name, $type)
{
	switch($type)
	{
		case FILE_DATASET_FAS :
			$dir = DIR_DATASET_FAS;
			$extension = "fas";
			break;
		case FILE_DATASET_QUAL :
			$dir = DIR_DATASET_QUAL;	
			$extension = "qual";
			break;
		case FILE_REFSET_FAS :
			$dir = DIR_REFSET_FAS;	
			$extension = "fas";
			break;
		case FILE_TAGSET_TXT :
			$dir = DIR_TAGSET_TXT;	
			$extension = "txt";
			break;
	        case FILE_FILE :
		        $dir = DIR_FILE;
			$extension = "dat";
			break;
			  
		default :
			$dir = DIR_OTHER;
			$extension = "dat";

	}
	
	return move_uploaded_file($_FILES[$name]["tmp_name"], "{$dir}/{$new_name}.{$extension}");
}

/**
 * Skickar vidare användare till sida
 * @param String Adress till sida
 * @param int Fördröjning, i sekunder. Används bara om headers är skickade.
 * @return void
 */
function redirect($url, $delay = 0)
{
	if (headers_sent())
	{
		echo "
		<script type=\"text/javascript\">
			function redirect_to(url)
			{
				window.location = '{$url}';
			}
		setTimeout(\"redirect_to('{$url}')\", {$delay}*1000);
		</script>";
	}
	else
	{
		header("location:{$url}");
	}
}

/**
 * Genererar slumpad, alfanumerisk sträng
 * @param int Längd, default 32
 * @return string Slumpad sträng
 */
function rand_string($lenth = 32) { 
    $aZ09 = array_merge(range('A', 'Z'), range('a', 'z'),range(0, 9)); 
    $out =''; 
    for($c=0; $c < $lenth; $c++)
	{ 
		$out .= $aZ09[mt_rand(0,count($aZ09)-1)]; 
    } 
    return $out; 
}

/**
 * Hämtar id för alla parameterset
 * @return int Array med ids
 */
function get_all_parametersets()
{
	$db = dbObject::get_instance();
	$user = new User();
	if ($user->is_admin())
	{
		$query = "SELECT idSettingsSet FROM ParameterSets";
	}
	else
	{
		$query = "SELECT idSettingsSet FROM ParameterSets WHERE owner = " . $user->get_id();
	}
	$db->query($query);
	$parametersets = Array();
	while ($parameterset = $db->fetch_row())
	{
		$parametersets[] = $parameterset[0];
	}

	return $parametersets;
}

/**
 * Hämtar id för alla parametrar
 * @return int Array med ids
 */
function get_all_parameters()
{
	$db = dbObject::get_instance();
	$query = "SELECT idParameter FROM Parameters ORDER BY parameter_order";
	$db->query($query);
	$parameters = Array();
	while ($parameter = $db->fetch_row())
	{
		$parameters[] = $parameter[0];
	}

	return $parameters;
}

/**
 * Hämtar id för alla tagsets
 * Endast tagsets som användare kan se
 */
function get_all_tagsets()
{
	$db = dbObject::get_instance();
	$user = new User();
	if ($user->is_admin())
	{
		$query = "SELECT idTagSets FROM TagSets";
	}
	else
	{
		$query = "SELECT idTagSets FROM TagSets WHERE owner = " . $user->get_id();
	}
	$db->query($query);
	$tagsets = Array();
	while ($tagset = $db->fetch_row())
	{
		$tagsets[] = $tagset[0];
	}

	return $tagsets;
}

/**
 * Hämtar id för alla files
 * Endast filer som användare kan se
 */
function get_all_files()
{
	$db = dbObject::get_instance();
	$user = new User();
	if ($user->is_admin())
	{
		$query = "SELECT idFiles FROM Files";
	}
	else
	{
		$query = "SELECT idFiles FROM Files WHERE owner = " . $user->get_id();
	}
	$db->query($query);
	$files = Array();
	while ($file = $db->fetch_row())
	{
		$files[] = $file[0];
	}

	return $files;
}


/**
 * Hämtar id för alla reference sequences
 * Endast reference sequences som användare kan se
 */
function get_all_refseqs()
{
	$db = dbObject::get_instance();
	$user = new User();
	if ($user->is_admin())
	{
		$query = "SELECT idReferenceSet FROM ReferenceSets";
	}
	else
	{
		$query = "SELECT idReferenceSet FROM ReferenceSets WHERE owner = " . $user->get_id() . " OR isPublic = 1";
	}
	$db->query($query);
	$refseqs = Array();
	while ($refseq = $db->fetch_row())
	{
		$refseqs[] = $refseq[0];
	}

	return $refseqs;
}

/**
 * Hämtar id för alla dataset
 * Endast dataset som användare kan se
 */
function get_all_datasets()
{
	$db = dbObject::get_instance();
	$user = new User();
	if ($user->is_admin())
	{
		$query = "SELECT idDatasets FROM Datasets";
	}
	else
	{
		$query = "SELECT idDatasets FROM Datasets WHERE owner = " . $user->get_id();
	}
	$db->query($query);
	$datasets = Array();
	while ($dataset = $db->fetch_row())
	{
		$datasets[] = $dataset[0];
	}

	return $datasets;
}

function get_all_datasets_ready()
{
        $db = dbObject::get_instance();
        $user = new User();
        if ($user->is_admin())
        {
                $query = "SELECT idDatasets FROM Datasets";
        }
        else
        {
                $query = "SELECT idDatasets FROM Datasets WHERE ready = 1 AND owner = " . $user->get_id();
        }
        $db->query($query);
        $datasets = Array();
        while ($dataset = $db->fetch_row())
        {
                $datasets[] = $dataset[0];
        }

        return $datasets;
}


/**
 * Hämtar id för alla job
 * Endast job som användare kan se
 */
function get_all_jobs()
{
	$db = dbObject::get_instance();
	$user = new User();
	if ($user->is_admin())
	{
		$query = "SELECT idJobs FROM Jobs";
	}
	else
	{
		$query = "SELECT idJobs FROM Jobs WHERE owner = " . $user->get_id();
	}
	$db->query($query);
	$jobs = Array();
	while ($job = $db->fetch_row())
	{
		$jobs[] = $job[0];
	}

	return $jobs;
}

/**
 * Meddelar användare genom något som påkallar uppmärksamhet
 * @param string Meddelande
 */
function notify($message)
{
	echo "
	<script type=\"text/javascript\">
		window.alert('Notification: {$message}');
	</script>";
}
?>
