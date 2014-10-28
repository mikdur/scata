<?php
if (isset($_GET['message']))
{
	$messages = Array();
	$messages['deleted'] = "File has been deleted";
	$messages['added'] = "File uploaded";
	notify($messages[$_GET['message']]);
}
?>

<h1>Upload</h1>
<span class="medium"><a href="./?p=upload&amp;do=add">Upload files</a></span><br />
<?php

echo "
<table>
	<thead>
		<tr>
			<th>Name</th><th>Size</th><th>Delete</th>
		</tr>
	</thead>
	<tbody>";
	
$files = get_all_files();
$form_num = 0;
if (count($files) > 0)
{
	foreach($files as $file)
	{
		$file = new File($file);
		echo "
			<tr>
				<td>{$file->get_name()}</td>
                                <td></td>
				<td>
					<form method=\"post\" action=\"./?p=upload&do=delete\" id=\"form_{$form_num}\">
						<input type=\"hidden\" value=\"{$file->id}\" name=\"id\">
						<a href=\"#\" onClick=\"confirm_delete({$form_num});\"><img src=\"./images/delete.png\" alt=\"delete\" /></a>
					</form>
				</td>
			</tr>";
		$form_num++;
	}
}
else
{
	echo "<tr><td>No files found</td></tr>";
}
echo "
	</tbody>
</table>";
?>
