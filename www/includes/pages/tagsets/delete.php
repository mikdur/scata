<h2>Delete</h2>
<?php
if (isset($_POST['id']))
{
	$tagset = new Tagset($_POST['id']);
	$name = $tagset->get_name();
	if ($tagset->delete())
	{
		redirect("./?p=tagsets&message=deleted");
	}
	else
	{
		echo "Couldn't delete tag set {$name}. Make sure you have the right permissions and that no job is using the set.<br />";
	}
}
else
{
	redirect("./?p=tagsets");
}

?>