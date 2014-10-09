<?php
if (isset($_GET['message']))
{
	$messages = Array();
	$messages['deleted'] = "Reference sequence has been deleted";
	$messages['added'] = "Added reference sequence";
	notify($messages[$_GET['message']]);
}
?>

<h1>Reference sequences</h1>
<span class="medium"><a href="./?p=refseqs&amp;do=add">Upload new set of reference sequences</a></span><br />
<?php

echo "
<table>
	<thead>
		<tr>
			<th>Name</th><th>Description</th><th>Status</th><th>Delete</th>
		</tr>
	</thead>
	<tbody>";
	
$refseqs = get_all_refseqs();
$form_num = 0;
foreach($refseqs as $refseq)
{
	$refseq = new Refseq($refseq);
	echo "
		<tr>
			<td>{$refseq->get_name()}</td>
                        <td>{$refseq->get_description()}</td>
			<td>".($refseq->is_ready() ? "OK" : "pending validation")."
			".($refseq->is_public() ? "Public sequence" : " ")."
			".($refseq->is_owner() ? " " : "provided by other user")."</td>
			<td>
				<form method=\"post\" action=\"./?p=refseqs&do=delete\" id=\"form_{$form_num}\">
					<input type=\"hidden\" value=\"{$refseq->id}\" name=\"id\">
					<a href=\"#\" onClick=\"confirm_delete({$form_num});\"><img src=\"./images/delete.png\" alt=\"delete\" /></a>
				</form>
			</td>
		</tr>";
	$form_num++;
}

echo "
	</tbody>
</table>";
?>
