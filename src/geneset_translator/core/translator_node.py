"""Translator graph node and edge dataclasses.

These classes represent nodes and edges from the NCATS Translator system,
including annotations from the Node Annotator API.
"""

from dataclasses import dataclass
from typing import Any, List, Optional


@dataclass
class TranslatorAttribute:
    """Translator node or edge attribute.

    Attributes are key-value pairs that provide additional metadata
    about nodes or edges in the Translator knowledge graph.
    """

    attribute_type_id: str
    """The type/key of this attribute (e.g., 'type_of_gene', 'taxid')."""

    value: Any
    """The value of this attribute."""

    value_type_id: Optional[str] = None
    """Optional type identifier for the value."""

    original_attribute_name: Optional[str] = None
    """Original name of the attribute from source database."""

    value_url: Optional[str] = None
    """Optional URL providing more information about the value."""

    attribute_source: Optional[str] = None
    """Source database or API that provided this attribute."""

    description: Optional[str] = None
    """Human-readable description of this attribute."""

    attributes: Optional[list] = None
    """Nested attributes (for complex attribute structures)."""


@dataclass
class TranslatorNode:
    """Translator graph node.

    Represents a node in the Translator knowledge graph, such as a gene,
    disease, protein, or biological process. Includes annotations from
    the Node Annotator API.
    """

    curie: str
    """CURIE identifier (e.g., 'NCBIGene:348', 'MONDO:0004975')."""

    label: Optional[str] = None
    """Human-readable name for the node."""

    types: Optional[List[str]] = None
    """List of biolink types (e.g., ['biolink:Gene', 'biolink:NamedThing'])."""

    synonyms: Optional[List[str]] = None
    """List of synonymous labels/names."""

    curie_synonyms: Optional[List[str]] = None
    """List of synonymous CURIE identifiers."""

    attributes: Optional[List[TranslatorAttribute]] = None
    """List of node attributes from Node Annotator."""

    taxa: Optional[List[str]] = None
    """List of taxa for this node (e.g., ['NCBITaxon:9606'] for human)."""

    @property
    def identifier(self) -> str:
        """Alias for curie - the node's identifier."""
        return self.curie

    @identifier.setter
    def identifier(self, value: str) -> None:
        """Set the node's identifier."""
        self.curie = value

    @property
    def categories(self) -> Optional[List[str]]:
        """Alias for types - the node's biolink categories."""
        return self.types

    @classmethod
    def from_dict(cls, data_dict: dict, return_synonyms: bool = False) -> "TranslatorNode":
        """Create a TranslatorNode from a dictionary.

        Args:
            data_dict: Dictionary containing node data. Must have 'curie' key.
            return_synonyms: If True, also extract synonyms from the dict.

        Returns:
            TranslatorNode instance.

        Raises:
            ValueError: If 'curie' key is missing from data_dict.
        """
        if 'curie' not in data_dict:
            raise ValueError('The input data dict must have a "curie" key.')

        node = cls(data_dict['curie'])

        if 'label' in data_dict:
            node.label = data_dict['label']

        if 'types' in data_dict:
            # Ensure types have biolink: prefix
            node.types = [
                f"biolink:{ty}" if not ty.startswith('biolink:') else ty
                for ty in data_dict['types']
            ]

        if 'taxa' in data_dict:
            node.taxa = data_dict['taxa']

        if return_synonyms:
            if 'synonyms' in data_dict:
                node.synonyms = data_dict['synonyms']
            elif 'names' in data_dict:
                # NameRes refers to synonyms as "names"
                node.synonyms = data_dict['names']

        return node

    def get_attribute(self, attribute_type_id: str) -> Optional[TranslatorAttribute]:
        """Get an attribute by its type ID.

        Args:
            attribute_type_id: The type/key of the attribute to find.

        Returns:
            The matching TranslatorAttribute, or None if not found.
        """
        if not self.attributes:
            return None
        for attr in self.attributes:
            if attr.attribute_type_id == attribute_type_id:
                return attr
        return None

    def get_attribute_value(self, attribute_type_id: str, default: Any = None) -> Any:
        """Get an attribute's value by its type ID.

        Args:
            attribute_type_id: The type/key of the attribute.
            default: Value to return if attribute not found.

        Returns:
            The attribute's value, or default if not found.
        """
        attr = self.get_attribute(attribute_type_id)
        return attr.value if attr else default


@dataclass
class TranslatorEdge:
    """Translator graph edge.

    Represents an edge/relationship between two nodes in the Translator
    knowledge graph.
    """

    subject: str
    """The subject node CURIE (source of the relationship)."""

    object: str
    """The object node CURIE (target of the relationship)."""

    predicate: str
    """The biolink predicate describing the relationship."""

    sources: Optional[list] = None
    """List of sources/databases that support this edge."""

    attributes: Optional[List[TranslatorAttribute]] = None
    """List of edge attributes."""
