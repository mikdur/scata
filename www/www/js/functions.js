function confirm_delete(form_id)
{
	var do_delete;
	do_delete = confirm("Are you sure you want to delete this item?");
	if (do_delete)
	{
		document.forms[form_id].submit();
	}
	else
	{
		return false;
	}
}