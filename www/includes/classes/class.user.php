<?php
/**
 * Klass för hantering av användare.
 * Sköter inloggning och redigering av uppgifter
 *
 */


class User
{
	private $db;
	public function __construct()
	{
		$this->db = dbObject::get_instance();
	}
	
	/**
	 * Logga in användare
	 * @param string E-postadress
	 * @param string Lösenord
	 * @return boolean True om inloggning lyckas, false annars.
	 */
	public function login($email, $password)
	{
	  if(preg_match("/[^@]+@[^@]+/", $email)){
	    $this->db->query("SELECT Password FROM Users WHERE Email = '{$email}'");
	    if ($this->db->row_count == 0)
	      {
		return USER_INVALID_CREDENTIALS;
	      }
	    
	    $pwd = $this->db->fetch_data();	
	    if (md5($password) != $pwd )
	      {
		if(md5($password . PASSWORD_SALT) != $pwd){
		     return USER_INVALID_CREDENTIALS;
		}else{
		  $npwd = md5($password);
		  $this->db->query("UPDATE Users SET Password = '{$npwd}' where Email = '{$email}'");
		}
	      }
	  }else{
	    /* AD user */
	       if($_SERVER['HTTPS'] != "on" && !pam_auth($email,$password)){
		 return USER_INVALID_CREDENTIALS;
	       }
	       
	       $email = $email . "@my-mgrid.mykopat.slu.se";
	       
	       $this->db->query("SELECT Password FROM Users WHERE Email = '{$email}'");
	       if ($this->db->row_count == 0)
		 {
		   $this->db->query("INSERT INTO Users (idUsers, Email, Password, Modified, Created, Enabled, isAdmin, sgeProject) VALUES(NULL, '${email}', NULL, NOW(), NOW(), 1, 0, 'mykopat-scata_slu')");
		 }
	  }
	  
	  $this->db->query("SELECT idUsers, Enabled, isAdmin FROM Users WHERE Email = '{$email}'");
	  $data = $this->db->fetch_assoc();
	  
	  if ($data['Enabled'] == 0)
	    {
	      return USER_NOT_ENABLED;
	    }
	  $_SESSION['uid'] = $data['idUsers'];
	  $_SESSION['admin'] = $data['isAdmin'];
	  
	  return USER_SUCCESS;
	}
	  
	/**
	 * Logga ut användare
	 */
	public function logout()
	{
		$_SESSION['uid'] = 0;
		$_SESSION['admin'] = 0;
		return USER_SUCCESS;
	}
	
	/**
	 * Se om användare är inloggad
	 * @return bool True om användare inloggad, annars false
	 */
	public function is_logged_in()
	{
		return (isset($_SESSION['uid']) && $_SESSION['uid'] != 0);
	}
	
	/**
	 * Se om användare har admin-rättigheter
	 * @return boolean True om den har admin-rättigheter, annars false.
	 */
	public function is_admin()
	{
		return (isset($_SESSION['admin']) && $_SESSION['admin'] == 1);
	}
	
	/**
	 * Hämta användares id
	 * @return int id
	 */
	public function get_id()
	{
		return $this->is_logged_in() ? $_SESSION['uid'] : 0;
	}
	
/** Funktioner för lösenordsbyte **/

	/**
	 * Begär byte av lösenord
	 * @param string Användares e-post
	 * @return boolean False om användaren inte finns, true annars.
	 */
	public function request_password_change($email)
	{
		// Verifiera att användare finns
		if ($email == "")
		{
			return 0;
		}
		$query = "SELECT COUNT(idUsers) FROM Users WHERE Email = '{$email}'";
		$this->db->query($query);
		if ($this->db->fetch_data() == 0)
		{
			return false;
		}
		
		$query = "SELECT
				idUsers
			FROM
				Users
			WHERE
				Email = '{$email}'";
		$this->db->query($query);
		$id = $this->db->fetch_data();

		// Ta bort gamla verifieringar
		$execute = "DELETE FROM UsersNewPassword WHERE user = {$id}";
		$this->db->execute($execute);

		// Skriv ny
		$verificationcode = md5($id . date("Y-m-d H:i") . PASSWORD_SALT);
		
		$insert = "INSERT INTO
			UsersNewPassword(
					user,
					verificationCode
			)
			VALUES
			(
				{$id},
				'{$verificationcode}'
			)";
			$this->db->execute($insert);
			
			// Skicka e-mail
			$subject = "SCATA: Password reset";
			$body = "Here's the link you need in order to reset your account password. Please note that it's only valid for 48 hours.\n\n".SITE_URL."/?p=lost_password&uid={$id}&verificationcode={$verificationcode}\n\nIf it's not clickable, copy and paste it in your browser's address field.\n\nIf you for some reason didn't request to reset your password, please ignore this e-mail. Your password has not yet been changed and everything is in order.";
			$body = wordwrap($body, 75);
			$headers = "From: SCATA <noreply@mykopat.slu.se>\r\nContent-type: text/plain; Charset=UTF-8\r\n";

			if (!mail($email, $subject, $body, $headers))
			{
				$this->delete($this->uid);
				throw new Exception("User: request_password_Change(): Couldn't send mail ({$email}).");
			}
			return true;
	}

	/**
	 * Se om byte för användare ska ske.
	 * Sant om inlägg i tabell finns och om skillnaden med dess timestamp är mindre än 48h.
	 * @param int Användarid
	 * @return boolean True om byte ska ske, false annars.
	 */
	public function is_password_to_be_changed($id)
	{
		$query = "SELECT
				(COUNT(idNewPassword) > 0) AND (HOUR(TIMEDIFF(requested, NOW())) < 48)
			FROM
				UsersNewPassword
			WHERE
				user = {$id}
			ORDER BY
				requested DESC
			LIMIT 1";
		$this->db->query($query);
		$result = $this->db->fetch_data();
		
		return ($result == "1");
	}
	
	/**
	 * Validerar verifieringskod
	 * @param int Användarid
	 * @param string Verifieringskod
	 * @return boolean True om koden är korrekt, false annars.
	 */
	public function validate_verification_code($id, $code)
	{
		if (!$this->is_password_to_be_changed($id))
		{
			throw new Exception("User: validate_verification_code(): Password not to be changed.");
		}
		$query = "SELECT
				COUNT(idNewPassword) = 1
			FROM
				UsersNewPassword
			WHERE
				user = {$id}";
		$this->db->query($query);
			
		return ($this->db->fetch_data() == 1);
	}
	
	/**
	 * Ger användaren nytt lösenord
	 * @param int Användarid
	 * @param string Lösenord i klartext
	 */
	public function change_password($id, $password)
	{
		if (!$this->is_password_to_be_changed($id))
		{
			throw new Exception("User: change_password(): Password not to be changed.");
		}
		$password = Registration::generate_password($password);
		$query = "UPDATE Users
			SET
				Password = '{$password}'
			WHERE
				idUsers = {$id}";
		$this->db->execute($query);
		$execute = "DELETE FROM
				UsersNewPassword
			WHERE user = {$id}";
		$this->db->execute($execute);
	}
}
?>
