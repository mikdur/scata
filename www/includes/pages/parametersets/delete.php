<h2>Delete</h2>
<?php
if (isset($_POST['id']))
{
	$paramset = new ParameterSet($_POST['id']);
	$name = $paramset->get_name();
	if ($paramset->delete())
	{
		redirect("./?p=parametersets&message=deleted");
	}
	else
	{
		echo "Couldn't delete parameter set {$name}. Make sure you have the right permissions and that no job is using the set.<br />";
	}
}
else
{
	redirect("./?p=parametersets");
}

?>