<?php
/**
 * Klass för registrering av användare
 */

class Registration
{
	private $db;
	private $uid;
	public $email;
	public $password;
	
	public function __construct($uid = 0, $email = "")
	{
		$this->db = dbObject::get_instance();
		$this->uid = $uid;
		$this->email = $email;
		
		if ($this->uid != 0)
		{
			//Hämta e-post
			$query = "SELECT Email, Password FROM Users WHERE idUsers = {$this->uid}";
			$this->db->query($query);
			if ($this->db->row_count == 0)
			{
				throw new Exception("Registration: __construct(): User id doesn't exist.");
			}
			$result = $this->db->fetch_assoc();
			$this->email = $result['Email'];
            $this->password = $result['Password'];
		}
		elseif ($email != "")
		{
			$query = "SELECT idUsers, Password FROM Users WHERE Email = '{$email}'";
			$this->db->query($query);
			if ($this->db->row_count == 0)
			{
				throw new Exception("Registration: __construct(): User id doesn't exist.");
			}
			$result = $this->db->fetch_assoc();
			$this->uid = $result['idUsers'];
            $this->password = $result['Password'];
		}
	}
	
	/**
	 * Verifierar att e-post inte redan finns i databasen
	 * Verifierar att e-post är giltig
	 * @param string Första angivna e-postadressen
	 * @param string Andra angivna e-postadressen
	 * @return int REGISTRATION_USER_EXISTS om användaren är registrerad.
	 */
	public function check_email($email_one, $email_two)
	{
		// Kontrollera att adresserna är lika
		if ($email_one != $email_two)
		{
			return REGISTRATION_NOT_MATCHING_EMAILS;
		}
		
		// Kontrollera att adress inte redan finns
		$this->db->query("SELECT COUNT(idUsers) FROM Users WHERE Email = '{$email_one}'");
		$count = $this->db->fetch_data();
		if ($count > 0)
		{
			return REGISTRATION_EMAIL_EXISTS;
		}
	
		
		// Annars är allt lugnt
		$this->email = $email_one;
		return REGISTRATION_SUCCESS;
	
	}
	
	/**
	 * Verifierar att lösenord är giltigt
	 * @param string Första angivna lösenordet
	 * @param string Andra angivna lösenordet
	 * @return array REGISTRATION_NOT_MATCHING_PASSWORDS om lösenordet inte matchar, REGISTRATION_NOT_VALID_PASSWORD om lösenordet är kortare än sex tecken. Annars REGISTRATION_SUCCESS.
 	 */
	public function check_passwords($password_one, $password_two)
	{
		// Kontrollera lika lösenord
		if ($password_one != $password_two)
		{
			return REGISTRATION_NOT_MATCHING_PASSWORDS;
		}
		
		// Kontrollera lösenord längre än sex tecken
		if (strlen($password_one) < 6)
		{
			return REGISTRATION_NOT_VALID_PASSWORD;
		}
		
		// Annars är allt lugnt
		$this->generate_password($password_one);
		return REGISTRATION_SUCCESS;
	}
	
	/**
	 * Krypterar lösenord
	 * @param string Lösenord i klartext
	 * @return string Krypterat lösenord
	 */
	public function generate_password($password)
	{
		$password = md5($password . PASSWORD_SALT);
		$this->password = $password;
		return $password;
	}
	
	/**
	 * Skapa användare i databasen
	 * @return void
	 */
	public function write_to_db()
	{
		if (!isset($this->email) || $this->email == "")
		{
			throw new Exception("Registration: write_to_db(): E-mail not set.");
		}
		if (!isset($this->password) || $this->password == "")
		{
			throw new Exception("Registration: write_to_db(): Password not set.");
		}

		$this->db->execute("INSERT INTO
			Users(
				Email,
				Password,
				Created,
				Enabled)
			VALUES(
				'{$this->email}',
				'{$this->password}',
				NOW(),
				false)");
		$this->uid = $this->db->last_id;
	}
	
	/**
	 * Generera aktiveringsnyckel
	 * @return string Aktiveringsnyckel
	 */
	function generate_activation_key()
	{
		if (!isset($this->uid) || $this->uid == 0)
		{
			throw new Exception("Registration: generate_activation_key(): User id not set.");
		}
		if (!isset($this->password))
		{
			$this->db->query("SELECT Password FROM Users WHERE idUsers = {$this->uid}");
			if ($this->db->row_count == 0)
			{
				throw new Exception("Registration: generate_activation_key(): User doesn't exist.");
			}
			else
			{
				$this->password = $this->db->fetch_data();
			}
		}

		$key = md5($this->email . $this->password . PASSWORD_SALT);
		return $key;
	}

	/**
	 * Skicka mail till användaren med aktiveringslänk
	 * @param int Användarid
	 * @return void
	 */
	public function send_activation_mail()
	{
		if (!isset($this->email) || $this->email == "")
		{
			$query = "SELECT Email FROM Users WHERE idUsers = {$this->uid}";
			if (!isset($this->uid) || $this->uid == 0)
			{
				throw new Exception("Registration: send_activation_mail(): User id not set.");
			}
			$this->db->query($query);
			if ($this->db->row_count == 1)
			{
				$this->email = $this->db->fetch_data();
			}
			else
			{
				throw new Exception("Registration: send_activation_email(): User doesn't exist.");
			}
		}
		
		$subject = "SCATA: Activation mail";
		$body = "Welcome to SCATA.\nThis message contains the activation link for your user account. Click the link below to activate the account.\n".SITE_URL."/?p=register&do=activate&email={$this->email}&key=".$this->generate_activation_key();
		$body = wordwrap($body, 75);
		$headers = "From: SCATA <noreply@mykopat.slu.se>\r\nContent-type: text/plain; Charset=UTF-8\r\n";
		
		if (!mail($this->email, $subject, $body, $headers))
		{
			$this->delete($this->uid);
			throw new Exception("Registration: send_activation_email(): Couldn't send activation mail ({$this->email}).");
		}
	}
	
	/**
	 * Aktivera användare
	 * @param string Användarid
	 * @param string Aktiveringsnyckel
	 * @return boolean true om användarnyckel stämmer och användare aktiverad, false annars.
	 */
	public function activate_user($key)
	{
		if ($this->generate_activation_key() == $key)
		{
			$this->db->execute("UPDATE Users SET Enabled = 1 WHERE idUsers = {$this->uid}");
			return true;
		}
		else
		{
			return false;
		}
		return false;
	}
	
	public function register_user($email, $password)
	{
		if (!is_array($email) || count($email) != 2)
		{
			throw new Exception("Registration: register_user(): E-mail must be array with two entries.");
		}
		if (count($password) != 2 || count($email) != 2)
		{
			throw new Exception("Registration: register_user(): Password must be array with two entries.");
		}
		if ($this->check_email($email[0], $email[1]) != REGISTRATION_SUCCESS)
		{
			return $this->check_email($email[0], $email[1]);
		}
		if ($this->check_passwords($password[0], $password[1]) != REGISTRATION_SUCCESS)
		{
			return $this->check_passwords($password[0], $password[1]);
		}
		$this->write_to_db();
		$this->send_activation_mail();
		return REGISTRATION_SUCCESS;
	}
	
	/**
	 * Raderar användare
	 */
	public function delete()
	{
		if (!isset($this->uid) || $this->uid == 0)
		{
			throw new Exception("Registration: delete(): User id not set.");
		}
		$query = "DELETE FROM Users WHERE idUsers = {$this->uid}";
		$this->db->execute($query);
	}
}

?>