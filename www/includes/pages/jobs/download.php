<h2>Download job</h2>
<?php

if (isset($_GET['id']))
{
	$job = new Job($_GET['id']);
	// Verifiera att användare har rätt till filen
	if ($job->user_has_rights())
	{
		$url = "./results/{$job->get_filename()}";
		echo "Automatically redirecting you to file. If not redirected, <a href=\"{$url}\">click here.</a>";
		redirect($url);
	}
	else
	{
		$url = "./?p=jobs";
		echo "An error occurred. Redirecting you to jobs page. If not redirected, <a href=\"{$url}\">click here</a>.";
		redirect($url);
	}

}
else
{
	$url = "./?p=jobs";
	echo "An error occurred. Redirecting you to jobs page. If not redirected, <a href=\"{$url}\">click here</a>.";
	redirect($url);
}
?>
