<h2>Delete</h2>
<?php
if (isset($_REQUEST['id']))
{
	$file = new File($_POST['id']);
	$name = $file->get_name();
	if ($file->delete())
	{
		redirect("./?p=upload&message=deleted");
	}
	else
	{
		echo "Couldn't delete file {$name}.<br />";
	}
}
else
{
	redirect("./?p=file");
}

?>