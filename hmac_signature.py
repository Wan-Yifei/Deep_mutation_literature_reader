import hmac
import hashlib


def generate_hmac_signature(obj, key):
    """
    Generate HASH message signature for object.

    :param obj: bytes or bytearry object needs a HASH authentication;
    :param key: key for generating signature;
    :return: hmac_signature: HASH signature for provided object with specified key.
    """
    if not isinstance(key, (bytes, bytearray)):
        key = str.encode(key, encoding="utf8")
    hmac_signature = hmac.new(key, obj, hashlib.sha224).hexdigest()
    return hmac_signature


def verify_hmac_signature(obj, key, signature_old):
    """
    Verify HASH message signature.

    :param obj: bytes project used for generating signature;
    :param key: key for generating signature;
    :param signature_old: expected signature;
    :return: verification: bool, True (verification pass) or False (verification failed).
    """
    signature_new = generate_hmac_signature(obj, key)
    signature_old = signature_old.strip()  # remove potential \n, \r and so on.
    if hmac.compare_digest(signature_new, signature_old):
        verification = True
    else:
        verification = False
        print("Integrity check failed!")
    return verification
