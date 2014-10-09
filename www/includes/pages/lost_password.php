<?php
$form = new FormObject("");
$form->add_field("password", "password","New password", "password", "mandatory=true;confirm=true;confirm_label=Confirm new password;minlength=6;label_after=At least 6 letters");
$form->add_field("verificationcode", "verificationcode", "", "hidden");
$form->add_field("uid", "uid", "", "hidden");

$form_reset = new FormObject("");
$form_reset->add_field("username", "username", "E-mail address", "text", "mandatory=true;confirm=true;confirm_label=Confirm e-mail address");

$user = new User();

echo "<h1>Reset password</h1>";

/**
 * Skriv ut formulär för lösenordsbyte
 * @param int användarid
 * @param string verifieringskod
 */
function print_password_form($id, $verificationcode)
{
	global $form;
	$form->set_value("verificationcode", $verificationcode);
	$form->set_value("uid", $id);
	echo $form->generate_javascript();
	echo $form->generate_form("./?p=lost_password", "Change password", false);
}

if (isset($_GET['verificationcode']) && isset($_GET['uid']))
{
	$error = "";
	$verificationcode = $_GET['verificationcode'];
	$uid = $_GET['uid'];
	// Ska byta göras?
	if ($user->is_password_to_be_changed($uid))
	{
		if ($user->validate_verification_code($uid, $verificationcode))
		{
			// Visa formulär för byta lösenord
			print_password_form($uid, $verificationcode);

		}
		else
		{
			$error = "Incorrect verification code.";
		}
	}
	else
	{
		$error = "Password reset not requested.";
	}
	
	if ($error != "")
	{
		echo "<h2>Error</h2>\n
		{$error}";
	}
}
else
{
	if ($form->got_post())
	{
		$error = "";
		$form->handle_input_data();
		
		if (!$form->validate_form())
		{
			$error = "An error occurred. Please try filling out the form again.";
			print_password_form($form->get_value("uid"), $form->get_value("verificationcode"));
		}
		else
		{
			$form->handle_input_data();
			if (!$user->is_password_to_be_changed($form->get_value("uid")) || !$user->validate_verification_code($form->get_value("uid"), $form->get_value("verificationcode")))
			{
				redirect("./?p=lost_password");
			}
			else
			{
				$user->change_password($form->get_value("uid"), $form->get_value("password"));
				echo "<h2>Password changed</h2>
				Redirecting you to <a href=\"./?p=login\">login page</a>...<br /><br />";
				redirect("./?p=login", 5);
			}
		}
		if ($error != "")
		{
			echo "<h2>Error</h2>\n
			{$error}";
		}
	}
	elseif ($form_reset->got_post())
	{
		$form_reset->handle_input_data();
		
		if ($user->request_password_change($form_reset->get_value("username")))
		{
			echo "<h2>Password reset request sent</h2>
				An e-mail with instructions how to reset your password has been sent.";
		}
		else
		{
			echo "<h2>Error</h2>\n
			Username doesn't exist in database.";
		}
	}
	else
	{
		echo "Use the form below to request a password change. An e-mail with a confirmation link will be sent to you. Click that link to access the page where you can set your new password.";
		echo $form_reset->generate_javascript();
		echo $form_reset->generate_form("./?p=lost_password", "Request password", false);
	}
}
?>