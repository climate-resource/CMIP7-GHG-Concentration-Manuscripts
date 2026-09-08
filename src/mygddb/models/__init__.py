"""
Definition of all (data) models

The pro of doing it this way is that all models are defined together.
The con is that you have both database and SDK-like models in the same place.
This makes for a less obviously clear separation between these two ideas.
It also makes importing the SDK-like models more expensive,
because all the database model initialisation stuff happens at the same time.
"""
