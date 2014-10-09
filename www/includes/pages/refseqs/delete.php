<h2>Delete</h2>
<?php
if (isset($_POST['id']))
{
	$refseq = new Refseq($_POST['id']);
	$name = $refseq->get_name();
	if ($refseq->delete())
	{
		redirect("./?p=refseqs&message=deleted");
	}
	else
	{
		echo "Couldn't delete reference sequence {$name}. Make sure you have the right permissions and that no job is using the set.<br />";
	}
}
else
{
	redirect("./?p=refseqs");
}

?>