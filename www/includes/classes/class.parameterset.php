<?php
/**
 * Klass för parameterset
 */
class ParameterSet
{
	public $id;
	private $db;
	private $user;
	private $owner;
	private $name;
	private $description;

	public function __construct($id = 0)
	{
		$this->db = dbObject::get_instance();
		$this->user = new User();
		$this->id = $id;
	}
	
	/**
	 * Skriver set till databasen
	 * @return void
	 */
	public function write_to_db()
	{
		if (!isset($this->owner))
		{
			throw new Exception("ParameterSet: write_to_db(): Owner not set.");
		}
		if (!isset($this->name))
		{
			throw new Exception("ParameterSet: write_to_db(): Name not set.");
		}
		if (!isset($this->description))
		{
			throw new Exception("ParameterSet: write_to_db(): Description not set.");
		}

		$query = "INSERT INTO
			ParameterSets(
				createTime,
				owner,
				Name,
				Description)
			VALUES(
				NOW(),
				{$this->owner},
				'{$this->name}',
				'{$this->description}')";
				
		$this->db->execute($query);
		$this->id = $this->db->last_id;
	}
	
	/**
	 * Sätter namn
	 * @param string Namn
	 */
	public function set_name($name)
	{
		$this->name = $name;
	}
	
	/**
	 * Sätter beskrivning
	 * @param string Description
	 */
	public function set_description($description)
	{
		$this->description = $description;
	}
	
	/**
	 * Sätter ägare
	 * @param int Ägares id
	 */
	public function set_owner($owner)
	{
		$this->owner = $owner;
	}
	
	/**
	 * Ser om användare är ägare
	 * @return boolean True om ägare, false annars.
	 */
	public function is_owner()
	{
		if ($this->id == 0)
		{
			throw new Exception("ParameterSet: is_owner(): id not set");
		}
		if (!isset($this->owner))
		{
			$query = "SELECT owner FROM ParameterSets WHERE idSettingsSet = {$this->id}";
			$this->db->query($query);
			$this->owner = $this->db->fetch_data();
		}
		return $this->owner == $this->user->get_id();
	}
	
	/**
	 * Raderar set
	 * @return boolean True om settet togs bort, false om settet inte kan tas bort p g a beroende(n), eller ägandeskap.
	 */
	public function delete()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("ParameterSet: delete(): Id not set.");
		}
		// Verifiera att det inte finns några jobb som använder settet
		$query = "SELECT COUNT(idJobs) FROM Jobs WHERE jobSettingsSet = {$this->id}";
		$this->db->query($query);
		if($this->db->fetch_data() > 0)
		{
			return false;
		}
		
		// Se om användare har rättigheter att radera
		$query = "SELECT ".($this->user->is_admin() ? "1" : "0")." = 1 OR owner = ".$this->user->get_id()." FROM ParameterSets WHERE idSettingsSet = {$this->id}";
		$this->db->query($query);
		if ($this->db->fetch_data() == 0)
		{
			return false;
		}

		$query = "DELETE FROM ParameterSettoValue WHERE parameterSetId = {$this->id}";
		$this->db->execute($query);
		
		$query = "DELETE FROM ParameterSets WHERE idSettingsSet = {$this->id}";
		$this->db->execute($query);
		
		return true;
	}
	
	/**
	 * Hämtar namn
	 * @return string Namn
	 */
	public function get_name()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("ParameterSet: get_name(): Id not set.");
		}
		if (!isset($this->name))
		{
			$query = "SELECT name FROM ParameterSets WHERE idSettingsSet = {$this->id}";
			$this->db->query($query);
			$this->name = $this->db->fetch_data();
		}
		return $this->name;
	}
	/**
	 * Hämtar beskrivning
	 * @return string beskrivning
	 */
	public function get_description()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("ParameterSet: get_description(): Id not set.");
		}
		if (!isset($this->description))
		{
			$query = "SELECT description FROM ParameterSets WHERE idSettingsSet = {$this->id}";
			$this->db->query($query);
			$this->description = $this->db->fetch_data();
		}
		return $this->description;
	}
	
	/**
	 * Ser om användare har rättighet till parametersettet
	 * @return boolean True om användare har det, false annars.
	 */
	public function user_has_rights()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("ParameterSet: has_rights_owner(): Id not set.");
		}
		return ($this->is_owner() || $this->user->is_admin());
	}
}

?>