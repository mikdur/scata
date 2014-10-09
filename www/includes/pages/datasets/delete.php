<h2>Delete</h2>
<?php
if (isset($_POST['id']))
{
	$dataset = new Dataset($_POST['id']);
	$name = $dataset->get_name();
	if ($dataset->delete())
	{
		redirect("./?p=datasets&message=deleted");
	}
	else
	{
		echo "Couldn't delete dataset {$name}. Make sure you have the right permissions and that no job is using the set.<br />";
	}
}
else
{
	redirect("./?p=datasets");
}

?>