<?php
/**
 * Klass för hantering av tagsets
 */

class Tagset
{
	private $name;
	public $id;
	private $db;
	private $user;
	
	public function __construct($id = 0)
	{
		$this->db = dbObject::get_instance();
		$this->user = new User();
		$this->id = $id;
	}
	
	/**
	 * Sätt namn
	 * @param name Namn
	 * @return void
	 */
	public function set_name($name)
	{
		$this->name = $name;
	}
	
	/**
	 * Hämta namn
	 * @return string namn
	 */
	public function get_name()
	{
		if ($this->id == 0)
		{
			throw new Exception("Tagset: get_name(): id not set");
		}
		if (!isset($this->name))
		{
			$query = "SELECT tagSetName FROM TagSets WHERE idTagSets = {$this->id}";
			$this->db->query($query);
			$this->name = $this->db->fetch_data();
		}
		return $this->name;
	}
	
	/**
	 * Ser om tagset är redo
	 * @return boolean True om redo, false annars.
	 */
	public function is_ready()
	{
		if ($this->id == 0)
		{
			throw new Exception("Tagset: is_ready(): id not set");
		}
		$query = "SELECT ready FROM TagSets WHERE idTagSets = {$this->id}";
		$this->db->query($query);
		return $this->db->fetch_data();
	}

	/**
	 * Ser om användare är ägare
	 * @return boolean True om ägare, false annars.
	 */
	public function is_owner()
	{
		if ($this->id == 0)
		{
			throw new Exception("Tagset: is_owner(): id not set");
		}
		if (!isset($this->owner))
		{
			$query = "SELECT owner FROM TagSets WHERE idTagSets = {$this->id}";
			$this->db->query($query);
			$this->owner = $this->db->fetch_data();
		}
		return $this->owner == $this->user->get_id();
	}
	
		
	/**
	 * Skapa tagset
	 * Sätter $this->id till tagsettets id
	 * @return boolean True om det lyckades, false annars.
	 */
	public function create()
	{
		if ($this->user->get_id() == 0)
		{
			throw new Exception("Tagset: create(): user id not set");
		}
		if ($this->name == "")
		{
			throw new Exception("Tagset: create(): Name not set");
		}
		$query = "INSERT INTO TagSets(
			createTime,
			owner,
			tagSetName)
		VALUES(
			NOW(),
			".$this->user->get_id().",
			'{$this->name}')";
		$this->db->execute($query);
		
		$this->id = $this->db->last_id;
		
		return true;
	}
	
	/**
	 * Ta bort tagset
	 * @param int id
	 * @return void
	 */
	public function delete($id = 0)
	{
		if ($id == 0)
		{
			$id = $this->id;
		}
		if ($id == 0)
		{
			throw new Exception("Tagset: delete(): Id not set");
		}
		
		// Verifiera att det inte finns några jobb som använder settet
		$query = "SELECT tagSet FROM Jobs WHERE tagSet = {$this->id}";
		$this->db->query($query);
		if($this->db->fetch_data() > 0)
		{
			return false;
		}

		// Se om användare har rättigheter att radera
		$query = "SELECT ".($this->user->is_admin() ? 1 : 0)." = 1 OR owner = ".$this->user->get_id()." FROM TagSets WHERE idTagSets = {$this->id}";
		$this->db->query($query);
		if ($this->db->fetch_data() == 0)
		{
			return false;
		}
		$query = "DELETE FROM TagSets WHERE idTagSets = {$this->id}";
		
		// Radera filen
		return $this->db->execute($query) && @unlink(DIR_TAGSET_TXT ."/". $this->id . ".txt");
	}
	
	/**
	 * Ser om användare har rättighet till settet
	 * @return boolean True om användare har det, false annars.
	 */
	public function user_has_rights()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Tagset: user_has_rights(): Id not set.");
		}
		return ($this->is_owner() || $this->user->is_admin());
	}
}
?>