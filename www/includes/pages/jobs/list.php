<?php
if (isset($_GET['message']))
{
	$messages = Array();
	$messages['deleted'] = "Job has been deleted";
	$messages['added'] = "Added job";
	notify($messages[$_GET['message']]);
}
?>

<h1>Job List</h1>
<span class="medium"><a href="./?p=jobs&amp;do=add">Submit a new job</a></span><br />
<?php

echo "
<table>
	<thead>
		<tr>
			<th>Name</th>
			<th>Description</th>
			<th>Parameter set</th>
			<th>Status</th>
			<th>Download result</th>
			<th>Delete</th>
		</tr>
	</thead>
	<tbody>";
	
$jobs = get_all_jobs();
$form_num = 0;
if (count($jobs) > 0)
{
	foreach($jobs as $job)
	{
		$job = new Job($job);
		$parameterset = new ParameterSet($job->get_parameterset());
		echo "
			<tr id=\"job-{$job->id}\">
				<td><a href=\"./?p=jobs&do=details&id={$job->id}\">".$job->get_name()."</a></td>
				<td>".$job->get_description()."</td>
				<td><a href=\"./?p=parametersets&do=details&id={$parameterset->id}\">".$parameterset->get_name()."</td>
				<td>".$job->get_status()."</td>
				<td>".($job->is_ready() ? "<a href=\"./?p=jobs&do=download&id={$job->id}\">Download</a>" : "Not available")."</td>
				<td>
					<form method=\"post\" action=\"./?p=jobs&do=delete\" id=\"form_{$form_num}\">
						<input type=\"hidden\" value=\"{$job->id}\" name=\"id\"><a href=\"#\" onClick=\"confirm_delete({$form_num});\"><img src=\"./images/delete.png\" alt=\"delete\" /></a>
					</form>
				</td>
			</tr>";
			$form_num++;
	}
}
else
{
	echo "<tr><td>No jobs added</td></tr>";
}
echo "
	</tbody>
</table>";
?>
