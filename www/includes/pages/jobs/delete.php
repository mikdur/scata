<h2>Delete</h2>
<?php
if (isset($_POST['id']))
{
	$job = new Job($_POST['id']);
	$name = $job->get_name();
	if ($job->delete())
	{
		redirect("./?p=jobs&message=deleted");
	}
	else
	{
		echo "Couldn't delete job {$name}. Make sure you have the right permissions.<br />";
	}
}
else
{
	redirect("./?p=jobs");
}

?>