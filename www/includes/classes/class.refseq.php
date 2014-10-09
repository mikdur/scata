<?php
/**
 * Klass för hantering av reference sequences
 */

class Refseq
{
	private $name;
	public $id;
	private $db;
	private $user;
	private $description;
	private $cutseq;
	private $cutlen;
	private $public;
	private $owner;
	
	public function __construct($id = 0)
	{
		$this->db = dbObject::get_instance();
		$this->user = new User();
		$this->id = $id;
	}
	
	/**
	 * Sätt namn
	 * @param string Namn
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
			throw new Exception("refseq: get_name(): id not set");
		}
		if (!isset($this->name))
		{
			$query = "SELECT Name FROM ReferenceSets WHERE idReferenceSet = {$this->id}";
			$this->db->query($query);
			$this->name = $this->db->fetch_data();
		}
		return $this->name;
	}
	
	/**
	 * Sätt beskrivning
	 * @param string Beskrivning
	 * @return void
	 */
	public function set_description($description)
	{
		$this->description = $description;
	}
	
	/**
	 * Hämta namn
	 * @return string namn
	 */
	public function get_description()
	{
		if ($this->id == 0)
		{
			throw new Exception("refseq: get_description(): id not set");
		}
		if (!isset($this->description))
		{
			$query = "SELECT Description FROM ReferenceSets WHERE idReferenceSet = {$this->id}";
			$this->db->query($query);
			$this->description = $this->db->fetch_data();
		}
		return $this->description;
	}

	/**
	 * Sättcutseq
	 * @param string cutseq
	 * @return void
	 */
	public function set_cutseq($cutseq)
	{
		$this->cutseq = $cutseq;
	}
	
	/**
	 * Hämta namn
	 * @return string namn
	 */
	public function get_cutseq()
	{
		if ($this->id == 0)
		{
			throw new Exception("refseq: get_cutseq(): id not set");
		}
		if (!isset($this->cutseq))
		{
			$query = "SELECT CutSequence FROM ReferenceSets WHERE idReferenceSet = {$this->id}";
			$this->db->query($query);
			$this->cutseq = $this->db->fetch_data();
		}
		return $this->cutseq;
	}



	/**
	 * Sättcutlen
	 * @param string cutlen
	 * @return void
	 */
	public function set_cutlen($cutlen)
	{
		$this->cutlen = $cutlen;
	}
	
	/**
	 * Hämta namn
	 * @return string namn
	 */
	public function get_cutlen()
	{
		if ($this->id == 0)
		{
			throw new Exception("refseq: get_cutlen(): id not set");
		}
		if (!isset($this->cutlen))
		{
			$query = "SELECT CutLength FROM ReferenceSets WHERE idReferenceSet = {$this->id}";
			$this->db->query($query);
			$this->cutlen = $this->db->fetch_data();
		}
		return $this->cutlen;
	}



	/**
	 * Sätt huruvida set ska vara publikt eller inte
	 * @param boolean True om publikt, false annars.
	 * @return void
	 */
	public function set_public($public = true)
	{
		$this->public = $public;
	}
	
	/**
	 * Ser om set är publikt eller inte.
	 * @return boolean True om publikt, false annars.
	 */
	public function is_public()
	{
		if ($this->id == 0)
		{
			throw new Exception("refseq: is_public(): id not set");
		}
		if (!isset($this->public))
		{
			$query = "SELECT isPublic FROM ReferenceSets WHERE idReferenceSet = {$this->id}";
			$this->db->query($query);
			$this->public = $this->db->fetch_data();
		}
		return $this->public;
	}

	/**
	 * Ser om refset är redo
	 * @return boolean True om redo, false annars.
	 */
	public function is_ready()
	{
		if ($this->id == 0)
		{
			throw new Exception("refseq: is_ready(): id not set");
		}
		$query = "SELECT ready FROM ReferenceSets WHERE idReferenceSet = {$this->id}";
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
			throw new Exception("refseq: is_owner(): id not set");
		}
		if (!isset($this->owner))
		{
			$query = "SELECT owner FROM ReferenceSets WHERE idReferenceSet = {$this->id}";
			$this->db->query($query);
			$this->owner = $this->db->fetch_data();
		}
		return $this->owner == $this->user->get_id();
	}

		
	/**
	 * Skapa refset
	 * Sätter $this->id till refsettets id
	 * @return boolean True om det lyckades, false annars.
	 */
	public function create()
	{
		if ($this->user->get_id() == 0)
		{
			throw new Exception("refseq: create(): user id not set");
		}
		if ($this->name == "")
		{
			throw new Exception("refseq: create(): Name not set");
		}
		if (!isset($this->description) || $this->description == "")
		{
			throw new Exception("refseq: create(): Description not set");
		}
		if (!isset($this->cutseq))
		{
			throw new Exception("refseq: create(): Cutseq not set");
		}
		if (!isset($this->cutlen))
		{
			throw new Exception("refseq: create(): Cutlen not set");
		}
		if (!isset($this->public))
		{
			throw new Exception("refseq: create(): Public not set");
		}
		if ($this->public != 1)
		{
			$this->public = 0;
		}
		$query = "INSERT INTO ReferenceSets(
			createTime,
			owner,
			Name,
			Description,
                        CutSequence,
                        CutLength,
			isPublic
			)
		VALUES(
			NOW(),
			".$this->user->get_id().",
			'{$this->name}',
			'{$this->description}',
                        '{$this->cutseq}',
                        '{$this->cutlen}',
			{$this->public})";
		$this->db->execute($query);
		
		$this->id = $this->db->last_id;
		
		return true;
	}
	
	/**
	 * Ta bort
	 * Raderar både fil och inlägg i databasen
	 * @param int id
	 * @return boolean True om set raderat, false annars.
	 */
	public function delete($id = 0)
	{
		if ($id == 0)
		{
			$id = $this->id;
		}
		if ($id == 0)
		{
			throw new Exception("refseq: delete(): Id not set");
		}
		// Verifiera att det inte finns några jobb som använder settet
		$query = "SELECT COUNT(idJobs) FROM Jobs WHERE jobSettingsSet = {$this->id}";
		$this->db->query($query);
		if($this->db->fetch_data() > 0)
		{
			return false;
		}

		// Se om användare har rättigheter att radera
		$query = "SELECT ".($this->user->is_admin() ? "1" : "0")." = 1 OR owner = ".$this->user->get_id()." FROM ReferenceSets WHERE idReferenceSet = {$this->id}";
		$this->db->query($query);
		if ($this->db->fetch_data() == 0)
		{
			echo "No right";
			return false;
		}
		$query = "DELETE FROM ReferenceSets WHERE idReferenceSet = {$this->id}";
		
		// Radera filen
		return $this->db->execute($query) && @unlink(DIR_REFSEQ_FAS ."/". $this->id . ".fas");
	}
	
	/**
	 * Ser om användare har rättighet till sekvensen
	 * @return boolean True om användare har det, false annars.
	 */
	public function user_has_rights()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Job: is_owner(): Id not set.");
		}
		return ($this->is_owner() || $this->user->is_admin() || $this->is_public());
	}
}
?>