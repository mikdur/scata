<?php
/**
 * Klass för parametrar
 */

class Parameter
{
	public $id;
	private $type;
	private $value;
	private $description;
	private $set_id;
	private $val_id;
	
	public function __construct($id = 0, $set_id = 0)
	{
		$this->db = dbObject::get_instance();
		$this->set_id = $set_id;
		$this->id = $id;
	}
	
	/**
	 * Sätter id för set
	 * @param int Settets id
	 */
	public function set_set_id($set_id)
	{
		$this->set_id = $set_id;
	}

	/**
	 * Skriver parametern till databasen
	 * @return True om 
	 */
	public function write_to_db()
	{
		if (!isset($this->keyword))
		{
			throw new Exception("Parameter: write_to_db(): Keyword not set.");
		}
		if (!isset($this->name))
		{
			throw new Exception("Parameter: write_to_db(): Name not set.");
		}
		if (!isset($this->description))
		{
			throw new Exception("Parameter: write_to_db(): Description not set.");
		}
		if (!isset($this->type))
		{
			throw new Exception("Parameter: write_to_db(): Type not set.");
		}
		$query = "INSERT INTO Parameters(
			createTime,
			keyword,
			Name,
			Description,
			type)
		VALUES(
			NOW(),
			{$this->keyword},
			{$this->name},
			{$this->description},
			$this->type)";
		$this->db->execute($query);
		$this->id = $this->db->last_id;
	}
	
	/**
	 * Sätter id för värde
	 * @param int ID
	 */
	public function set_value_id($val_id)
	{
		$this->val_id = $val_id;
	}
	
	/**
	 * Sätter fördefinierat textvärde
	 * @param string Värde
	 * @return boolean True om värdet satt, false annars
	 */
	public function set_text_value($value)
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Parameter: set_text_value(): Id not set.");
		}
		$query = "INSERT INTO ParameterSettoValue(
			parameterSetId,
			parameterValueId,
			textValue)
		VALUES(
			{$this->set_id},
			".$this->get_default_value()->id.",
			'{$value}'
		)";
		$this->db->execute($query);
		$this->val_id = $this->db->last_id;
		return true;
	}
	
	/**
	 * Sätter selectvärde
	 * @param int Id för värde
	 * @return void
	 */
	public function set_select_value($value)
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Parameter: set_select_value(): Id not set.");
		}
		$query = "INSERT INTO ParameterSettoValue(
			parameterSetId,
			parameterValueId)
		VALUES(
			{$this->set_id},
			{$value}
				)";
		$this->db->execute($query);

	}
	
	/**
	 * Hämtar textvärde
	 * @param int ID för ParameterSet
	 * @return string värde
	 */
	public function get_text_value($set_id = 0)
	{
		if (isset($set_id) && $set_id != 0)
		{
			$query = "SELECT
				psv.textValue
			FROM
				ParameterSettoValue AS psv
				INNER JOIN
					ParameterValue AS pv
					ON
					pv.idValue = psv.parameterValueId
			WHERE
				psv.parameterSetId = {$set_id}
				AND
				pv.parameterId = {$this->id}";
		}

		elseif (isset($this->val_id) && $this->val_id != 0)
		{
			$query = "SELECT textValue FROM ParameterSettoValue WHERE parameterValueId = {$this->val_id}";
		}
		else
		{
			throw new Exception("Parameter: get_text_value(): Must either specify val-id or set-id");
		}
		$this->db->query($query);
		if ($this->db->row_count == 0)
		{
			throw new Exception("Parameter: get_text_value(): Id doesn't exist in db.");
		}

		return $this->db->fetch_data();
	}

	/**
	 * Hämtar select-värde
	 * @param int ID för ParameterSet
	 * @return ParameterValue Objekt som representerar värde
	 */
	public function get_select_value($set_id = 0)
	{
		if (isset($set_id) && $set_id != 0)
		{
			$query = "SELECT
				idValue
			FROM
				ParameterValue AS pv
				INNER JOIN
					ParameterSettoValue AS psv
					ON
						psv.parameterSetId = {$set_id}
			WHERE
				pv.idValue = psv.parameterValueId
				AND
				pv.parameterId = {$this->id}";
		}
		elseif (isset($this->val_id) && $this->val_id != 0)
		{
			$query = "SELECT parameterValueId FROM ParameterSettoValue WHERE parameterValueId = {$this->val_id}";
		}
		else
		{
			throw new Exception("Parameter: get_select_value(): Must either specify val-id or set-id");
		}
		$this->db->query($query);

		if ($this->db->row_count == 0)
		{
			throw new Exception("Parameter: get_select_value(): Id doesn't exist in db.");
		}

		return new ParameterValue($this->db->fetch_data(), $this->id);
	}

	/**
	 * Hämtar beskrivning
	 * @return string Beskrivning
	 */
	public function get_description()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Parameter: get_description(): Id not set.");
		}
		$query = "SELECT Description FROM Parameters WHERE idParameter = {$this->id}";
		$this->db->query($query);

		if ($this->db->row_count == 0)
		{
			throw new Exception("Parameter: get_description(): Id doesn't exist in db.");
		}

		return $this->db->fetch_data();
	}
	
	/**
	 * Sätter beskrivning
	 * @param string Beskrivning
	 * @return void
	 */
	public function set_description($description)
	{
		$this->description = $description;
	}

	/**
	 * Hämtar typ
	 * @return String "text" om typen är "text", "select" om typen är "select"
	 */
	public function get_type()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Parameter: get_type(): Id not set.");
		}
		$query = "SELECT type FROM Parameters WHERE idParameter = {$this->id}";
		$this->db->query($query);
		$this->type = $this->db->fetch_data();

		if ($this->db->row_count == 0)
		{
			throw new Exception("Parameter: get_type(): Id doesn't exist in db.");
		}

		return $this->type;
	}
	
	/**
	 * Sätter typ
	 * @param string typ
	 * @return boolean True om typen är satt till "text" eller "select". False annars.
	 */
	public function set_type($type)
	{
		if ($type != "text" || $type != "select")
		{
			throw new Exception("Parameter: set_type(): Incorrect type. Must be either \"select\" or \"text\".");
		}
		$this->type = $type;
		return true;
	}
	
	/**
	 * Hämtar samtliga möjliga värden för en parameter av typen "select"
	 * @return ParameterValue Array med värden
	 */
	public function get_possible_select_values()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Parameter: get_possible_select_values(): Id not set.");
		}

		$query = "SELECT idValue FROM ParameterValue WHERE parameterId = {$this->id}";


		$this->db->query($query);

		if ($this->db->row_count == 0)
		{
			throw new Exception("Parameter: get_possible_select_values(): Id doesn't exist in database.");
		}

		$parameter_values = Array();

		while ($param_id = $this->db->fetch_row())
		{
			$parameter_values[] = new ParameterValue($param_id[0], $this->id);
		}
		return $parameter_values;
		
	}
	
	/**
	 * Hämtar namn på parametern
	 * @return string Namn
	 */
	public function get_name()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Parameter: get_name(): Id not set.");
		}
		$query = "SELECT Name FROM Parameters WHERE idParameter = {$this->id}";
		$this->db->query($query);
		
		
		return $this->db->fetch_data();
	}
	
	/**
	 * Sätter namn
	 * @param string Namn
	 * @return void
	 */
	public function set_name($name)
	{
		$this->name = $name;
	}
	
	/**
	 * Ser om parametern är tillgänglig
	 * @return boolean True om tillgänglig, false annars
	 */
	public function is_available()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Parameter: get_name(): Id not set.");
		}
		$query = "SELECT available FROM Parameters WHERE idParameter = {$this->id}";
		$this->db->query($query);
		$retval = $this->db->fetch_data();
	
		return $retval == 1;
	}
	
	/**
	 * Hämta default-värde för parametern
	 * Finns inget default-värde returneras ett värde som är 0.
	 * @return ParameterValue Default-värde
	 */
	public function get_default_value()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("Parameter: get_default_value(): Id not set.");
		}
		$query = "SELECT idValue FROM ParameterValue WHERE isDefault = 1 AND parameterId = {$this->id}";
		$this->db->query($query);
		if ($this->db->row_count == 0)
		{
			$param_value = new ParameterValue(0, $this->id);
			$param_value->set_value("");
			$param_value->set_display_value("");
		}
		
		$id = $this->db->fetch_data();

		
		return new ParameterValue($id, $this->id);
	}
}
?>
