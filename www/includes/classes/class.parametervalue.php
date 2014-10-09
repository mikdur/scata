<?php

class ParameterValue
{
	public $id;
	public $param_id;
	private $db;
	private $value;
	private $display_value;

	public function __construct($id = 0, $param_id = 0)
	{
		$this->db = dbObject::get_instance();
		$this->id = $id;
		$this->param_id = $param_id;
		$this->display_value = "";
	}
	
	/**
	 * Hämtar displayvärde för parametervärde
	 * @return string Displayvärde
	 */
	public function get_display_value()
	{
		if (!isset($this->display_value) || $this->display_value == "")
		{
			if (!isset($this->id) || $this->id == 0)
			{
				throw new Exception("ParameterValue: get_display_value(): Id not set.");
			}
			$query = "SELECT displayValue FROM ParameterValue WHERE idValue = {$this->id}";
			$this->db->query($query);
			
			if ($this->db->row_count == 0)
			{
				throw new Exception("ParameterValue: get_display_value(): Id doesn't exist in database.");
			}
			$this->display_value = $this->db->fetch_data();
		}
		return $this->display_value;
	}
	
	/**
	 * Sätter displayvärde
	 * @param string Displayvärde
	 * @return void
	 */
	public function set_display_value($display_value)
	{
		$this->display_value = $display_value;
	}
	
	/**
	 * Sätter värde
	 * @param string Värde
	 */
	public function set_value($value)
	{
		$this->value = $value;
	}
	
	/**
	 * Hämtar värde
	 * @return string Värde
	 */
	public function get_value()
	{
		if (!isset($this->value))
		{
			if (!isset($this->id) || $this->id == 0)
			{
				$this->value = "";
			}
			else
			{
				$query = "SELECT value FROM ParameterValue WHERE idValue = {$this->id}";
				$this->db->query($query);

				if ($this->db->row_count == 0)
				{
					throw new Exception("ParameterValue: get_value(): Id doesn't exist in database.");
				}

				$this->value = $this->db->fetch_data();				
			}
		}
		return $this->value;
	}
	
	/**
	 * Sätter parametervärdet till default
	 * Kan endast köras om id är angivet
	 * @return void
	 */
	public function set_default()
	{
		if (!isset($this->param_id) || $this->param_id == 0)
		{
			throw new Exception("ParameterValue: set_default(): Parameter id not set.");
		}
		
		$query = "UPDATE ParameterValue SET isDefault = FALSE WHERE parameterId = {$this->param_id}";
		$this->db->execute($query);

		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("ParameterValue: set_default(): Id not set.");
		}

		$query = "UPDATE ParameterValue SET isDefault = TRUE WHERE idValue = {$this->id}";
		$this->db->execute($query);
	}
	
	/**
	 * Ser om parametervärdet är default
	 * Kan endast köras om id är angivet
	 * @return boolean True om parametervärdet är default, false annars
	 */
	public function is_default()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("ParameterValue: is_default(): Id not set.");
		}

		$query = "SELECT isDefault IS TRUE FROM ParameterValue WHERE idValue = {$this->id}";
		$this->db->execute($query);
		return $this->db->fetch_data() == 1;
	}
	
	/**
	 * Skapar nytt värde
	 * @return int ID för värdet
	 */
	public function create()
	{
		if (!isset($this->param_id) || $this->param_id == 0)
		{
			throw new Exception("ParameterValue: create(): Parameter id not set.");
		}
		if (!isset($this->value) || $this->value == "")
		{
			throw new Exception("ParameterValue: create(): Value not set.");
		}
		if (!isset($this->display_value) || $this->display_value == "")
		{
			throw new Exception("ParameterValue: create(): Parameter id not set.");
		}
		$query = "INSERT INTO ParameterValue(parameterId, value, displayValue) VALUES({$this->param_id}, '{$this->value}', '{$this->display_value}')";
		$this->db->execute($query);
		$this->id = $this->db->last_id;
		return $this->id;
	}

	/**
	 * Uppdaterar värde
	 */
	public function update()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("ParameterValue: update(): Id not set.");
		}
		if (!isset($this->value) || $this->value == "")
		{
			throw new Exception("ParameterValue: update(): Value not set.");
		}
		if (!isset($this->display_value) || $this->display_value == "")
		{
			throw new Exception("ParameterValue: update(): Parameter id not set.");
		}
		$query = "UPDATE ParameterValue
			SET
				value = {$this->value},
				displayValue = {$this->display_value}
			WHERE
				idValue = {$this->id}";
		$db->execute($query);
	}
	
	/**
	 * Raderar värde
	 */
	public function delete()
	{
		if (!isset($this->id) || $this->id == 0)
		{
			throw new Exception("ParameterValue: delete(): Id not set.");
		}
		$query = "DELETE FROM ParameterValue WHERE idValue = {$this->id}";
		$this->db->execute($query);
	}
}
?>